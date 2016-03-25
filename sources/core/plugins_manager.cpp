/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2014, 2015 CNRS
   Copyright (C) 2013, 2014, 2015, 2016 Inria

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
#include <cstdlib>
#include <stdio.h>
#include <list>

using namespace alta;

//#define DEBUG


//! Add dynamic library extension (.so or .dll) to a dynamic object as well as
//! prepending 'lib' depending on the plateform.
static std::string library_name(const std::string name)
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

// Return the plugin search path based on the 'ALTA_PLUGIN_PATH' environment
// variable.  Empty entries in that variable's value are taken to designate
// the plugin installation directory.
static std::list<std::string> plugin_search_path()
{
  std::list<std::string> dirs;
  std::string obj_str ;
  const char *env_str = std::getenv("ALTA_PLUGIN_PATH");

  if(env_str == NULL) {
    obj_str = "";
  } else {
    obj_str = env_str;
  }

  if (obj_str.empty())
      // Always add at least the system's plugin directory.
      dirs.push_back(ALTA_PLUGIN_DIRECTORY);
  else {
      for (auto start = 0, next = 0;
           start < obj_str.size();
           start = next + 1)
      {
          next = obj_str.find(':', start);
          if (next == std::string::npos)
              next = obj_str.size();

          // An empty entry denotes the system directory.
          auto entry = next == start
              ? ALTA_PLUGIN_DIRECTORY
              : obj_str.substr(start, next - start);

          dirs.push_back(entry);

          if (next == obj_str.size() - 1)
              // This is a trailing empty entry.
              dirs.push_back(ALTA_PLUGIN_DIRECTORY);
      }
  }

  return dirs;
}

//! \brief Open a dynamic library file (.so or .dll) and extract the associated
//! provide function. The template argument is used to cast the library to a
//! specific type.
template<typename T>
static T open_library(const std::string& filename, const char* function)
{
  auto directories = plugin_search_path();

  for (auto&& directory: directories)
  {
    auto libname = directory + "/" + library_name(filename);

#ifdef _WIN32
    HINSTANCE handle = LoadLibraryA(libname.c_str());
    if(handle != NULL)
    {
        T res = (T)GetProcAddress(handle, function);

        if(res == NULL)
        {
#ifdef DEBUG_CORE
        std::cerr << "<<ERROR>> unable to load the symbol \"" << function << "\" from " << filename << std::endl;
#endif
        continue;
        }
#ifdef DEBUG_CORE
        std::cout << "<<DEBUG>> will provide a " << function << " for library \"" << filename << "\"" << std::endl;
#endif
        return res;
    }
    else
    {
#ifdef DEBUG_CORE
        std::cerr << "<<ERROR>> unable to load the dynamic library file \"" << libname << "\"" << std::endl;
        std::cerr << "          cause: \"" << GetLastError() << "\"" << std::endl;
#endif
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
#ifdef DEBUG_CORE
            std::cerr << "<<ERROR>> unable to load the symbol \"" << function << "\" from " << libname << std::endl;
#endif
            continue;
        }
#ifdef DEBUG_CORE
        std::cout << "<<DEBUG>> will provide a " << function << " for library \"" << libname << "\"" << std::endl;
#endif
        return (T)res;
    }
    else
    {
#ifdef DEBUG_CORE
        std::cerr << "<<ERROR>> unable to load the dynamic library file \"" << libname << "\"" << std::endl;
      std::cerr << "          cause: \"" << dlerror() << "\"" << std::endl;
#endif
        continue;
    }
#endif
  }

  std::cerr << "<<ERROR>> unable to load the symbol \"" << function << "\" from " << filename << std::endl;
  return NULL;
}

//! \brief load a function from the ALTA input file.
function* plugins_manager::load_function(const std::string& filename)
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

    //RP: Returning NULL immediately now
    return NULL;
  }

  // Parse the header for the function command line and the dimension
  // of the function
  arguments header = arguments::parse_header(file);

  std::pair<int, int> dim = header.get_pair<int>("DIM");
  nX = dim.first;
  nY = dim.second;

  param_in = params::parse_input(header.get_string("PARAM_IN", "UNKNOWN_INPUT"));
  param_out = params::parse_output(header.get_string("PARAM_OUT", "UNKNOWN_OUTPUT"));
  args = arguments::create_arguments(header["CMD"]);

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
  if( f->load(file) )
  {
    return f;
  }

  std::cout << "<<ERROR>> COULD NOT LOAD THE BRDF from File " << filename << std::endl;
   return NULL;
}

//! Get an instance of the function selected based on the name <em>n</em>.
//! Return NULL if none exists.
function* plugins_manager::get_function(const arguments& args)
{
  #ifdef DEBUG
  std::cout << __FILE__ << " " << __LINE__ << std::endl;
  std::cout << " INitial args = " << args << std::endl;
  #endif

  if(!args.is_defined("func"))
  {
    #ifdef DEBUG
    std::cout << "<<DEBUG>> no function plugin specified, returning a rational function" << std::endl;
    #endif
    return new rational_function();
  }

  // The function to be returned.
  function* func = NULL;

  if(args.is_vec("func")) //Checking if  the argument prodived with --func is of the form --func [ arg1, arg2, arg3]
  {
    #ifdef DEBUG
    std::cout << __FILE__ << " " << __LINE__ << " ARGUMENTS of --func are forming a VECTOR " << std::endl;
    #endif 

    std::vector<std::string> args_vec = args.get_vec("func");

    // Treating the case []
    if(args_vec.size() == 0)
    {
        return NULL;
    }

   //! create a *compound* class to store multiple ! functions in it.
    compound_function* compound = new compound_function();

   //! For each args_vec element, create a function object and add ! it to the
   //compound one.
    for(unsigned int i=0; i<args_vec.size(); ++i)
    {
      std::string n("--func ");
      n.append(args_vec[i]);

      #ifdef DEBUG
      std::cout << __FILE__ << " " << __LINE__ << std::endl
                << " at i=" << i << " n=" << n << std::endl;
      std::cout << "<<DEBUG>> load function with args: " << n << std::endl;
      #endif
      
      arguments temp_args = arguments::create_arguments(n);

      //Recursive call
      function* f = get_function(temp_args);
      nonlinear_function *nl_f = dynamic_cast<nonlinear_function*>(f);
      if(nl_f == NULL)
      {
          std::cerr << "<<ERROR>> only non-linear functions are compatible with a compound" << std::endl;
      }
      else
      {
          compound->push_back(ptr<nonlinear_function>(nl_f), temp_args);
      }
    }

      //! return the compound class
      func = compound;
      #ifdef DEBUG
      std::cout << __FILE__ << " " << __LINE__ << " WE HAVE A COMPOUND " << std::endl;
      #endif
  }
  else //Here we check if the argument provided with --func is file describing a function
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
      #ifdef DEBUG
      std::cout << __FILE__ << " " << __LINE__ << std::endl;
      std::cout << " INSIDE fresnel args = " << args << std::endl;
      #endif

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

      #ifdef DEBUG
      std::cout << __FILE__ << " " << __LINE__ << std::endl;
      std::cout << " n = " << n << std::endl;
      #endif
       
      nonlinear_function* func_fres = dynamic_cast<nonlinear_function*>(get_function(arguments::create_arguments(n)));
      if(func_fres != NULL)
      {
        bool const fresnel_is_fixed = (args.is_defined("fixed-fresnel")) ? (true) : (false);
        bool const lobe_is_fixed    = (args.is_defined("fixed-lobe")) ? (true) : (false);
       
        if( fresnel_is_fixed )
        {
          std::cout << "<<DEBUG>> The Fresnel term is fixed" << std::endl;
        }
        if( lobe_is_fixed )
        {
          std::cout << "<<DEBUG>> The lobe is fixed" << std::endl;
        }

        return new product_function(ptr<nonlinear_function>(nl_func),
                                    ptr<nonlinear_function>(func_fres),
                                    lobe_is_fixed, fresnel_is_fixed);
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


ptr<function> plugins_manager::get_function(const std::string& n)
{
    if(n.empty())
    {
#ifdef DEBUG
        std::cout << "<<DEBUG>> no function plugin specified, returning nothing" << std::endl;
#endif
        return NULL;
    }

   FunctionPrototype myFunc = open_library<FunctionPrototype>(n, "provide_function");
    if(myFunc != NULL)
    {
#ifdef DEBUG
        std::cout << "<<DEBUG>> using function provider in file \"" << n << "\"" << std::endl;
#endif
        return ptr<function>(myFunc());
    }
    else
    {
        std::cerr << "<<ERROR>> no function provider found in file \"" << n << "\"" << std::endl;
        return NULL;
    }
}


ptr<data> plugins_manager::get_data(const std::string& n, const arguments& args)
{
    if(n.empty() || n == "vertical_segment")
    {
#ifdef DEBUG
        std::cout << "<<DEBUG>> no data plugin specified, returning a vertical_segment loader" << std::endl;
#endif
        return ptr<data>(new vertical_segment());
    }

   DataPrototype myData = open_library<DataPrototype>(n, "provide_data");
    if(myData != NULL)
    {
#ifdef DEBUG
        std::cout << "<<DEBUG>> using data provider in file \"" << n << "\"" << std::endl;
#endif
        return ptr<data>(myData(args));
    }
    else
    {
        std::cerr << "<<ERROR>> no data provider found in file \"" << n << "\"" << std::endl;
        return ptr<data>(new vertical_segment()) ;
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
        std::cout << "<<DEBUG>> using fitter provider in file \"" << n << "\"" << std::endl;
#endif
        return ptr<fitter>(myFitter());
    }
    else
    {
        std::cerr << "<<ERROR>> no fitter provider found in file \"" << n << "\"" << std::endl;
        return NULL ;
    }
}

void plugins_manager::check_compatibility( ptr<data>& d, 
                                           const ptr<function>& f,
                                           const arguments& args)
{
  if(d->parametrization().input_parametrization() == params::UNKNOWN_INPUT &&
    f->input_parametrization() == params::UNKNOWN_INPUT)
  {
    std::cout << "<<WARNING>> both function and data objects have no parametrization" << std::endl;
  }
  else
  {
    if(d->parametrization().input_parametrization() == params::UNKNOWN_INPUT)
    {
      std::cout << "<<WARNING>> unknown parametrization for data" << std::endl;
    }

    if(f->input_parametrization() == params::UNKNOWN_INPUT)
    {
      std::cout << "<<DEBUG>> function will take the parametrization of the data" << std::endl;
      f->setParametrization(d->parametrization().input_parametrization());
    }
    else if(d->parametrization().input_parametrization() != f->input_parametrization() && args.is_defined("change-param"))
    {
      std::cout << "<<INFO>> has to change the parametrization of the input data " << params::get_name(d->parametrization().input_parametrization()) << std::endl;
      std::cout << "<<INFO>> to " << params::get_name(f->input_parametrization()) << std::endl;
      ptr<data_params> dd = ptr<data_params>(new data_params(d, f->input_parametrization()));
      d = dynamic_pointer_cast<data>(dd) ;
    }
    else
    {
      std::cout << "<<DEBUG>> no change was made to the parametrization" << std::endl;
    }
  }

  if(f->dimY() != d->parametrization().dimY())
  {
    f->setDimY(d->parametrization().dimY());
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
// \todo : Why not move these functions to common.h ?
#ifdef _WIN32
#include <windows.h>
size_t plugins_manager::get_system_memory()
{
  MEMORYSTATUSEX status;
  status.dwLength = sizeof(status);
  GlobalMemoryStatusEx(&status);
  return status.ullTotalPhys;
}
#elif defined(__APPLE__)
#include <unistd.h>
#include <sys/types.h>
#include <sys/param.h>
#include <sys/sysctl.h>
size_t plugins_manager::get_system_memory()
{
  int mib[2];
  mib[0] = CTL_HW;
  mib[1] = HW_MEMSIZE;    
  int64_t size = 0;   
  size_t len = sizeof( size );
  if ( sysctl( mib, 2, &size, &len, NULL, 0 ) == 0 )
    return (size_t)size;
  
  return 0L;
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
