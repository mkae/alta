#include "plugins_manager.h"
#include "rational_function.h"
#include "vertical_segment.h"

/*
#include <QCoreApplication>
#include <QDir>
*/

#ifdef WIN32
    #include <windows.h>
#else
    #include <dlfcn.h>
#endif
#include <stdio.h>

//! \brief Open a dynamic library file (.so or .dll) and extract the associated
//! provide function. The template argument is used to cast the library to a
//! specific type.
template<typename T> T open_library(const std::string& filename, const char* function)
{
#ifdef WIN32
    HINSTANCE handle = LoadLibrary((LPCTSTR)filename.c_str());
    if(handle != NULL)
    {
        T res = (T)GetProcAddress(handle, function);

        if(res == NULL)
        {
            std::cerr << "<<ERROR>> unable to load the symbol \"" << function << "\" from " << filename << std::endl;
            return NULL;
        }
#ifdef DEBUG_CORE
        std::cout << "<<DEBUG>> will provide a " << function << " for library \"" << filename << "\"" << std::endl;
#endif
        return res;
    }
    else
    {
        std::cerr << "<<ERROR>> unable to load the dynamic library file \"" << filename << "\"" << std::endl;
        std::cerr << "          cause: \"" << GetLastError() << "\"" << std::endl;
        return NULL;
    }
#else
    void* handle = dlopen(filename.c_str(), RTLD_LAZY);
    if(handle != NULL)
    {
        void* res = dlsym(handle, function);

        if(dlerror() != NULL)
        {
            std::cerr << "<<ERROR>> unable to load the symbol \"" << function << "\" from " << filename << std::endl;
            return NULL;
        }
#ifdef DEBUG_CORE
        std::cout << "<<DEBUG>> will provide a " << function << " for library \"" << filename << "\"" << std::endl;
#endif
        return (T)res;
    }
    else
    {
        std::cerr << "<<ERROR>> unable to load the dynamic library file \"" << filename << "\"" << std::endl;
		  std::cerr << "          cause: \"" << dlerror() << "\"" << std::endl;
        return NULL;
    }
#endif
}


// Create the object, parse the argument and load all the plugins
plugins_manager::plugins_manager(const arguments& args) 
{
/*
	QDir pluginsDir;
	if(args.is_defined("plugins"))
	{
		pluginsDir = QDir(args["plugins"].c_str()) ;
	}
	else if(QCoreApplication::instance() != NULL)
	{
		pluginsDir = QDir(QCoreApplication::instance()->applicationDirPath());
	}
	*/

	/*
	foreach (QString fileName, pluginsDir.entryList(QDir::Files)) 
	{
		QPluginLoader loader ;
		loader.setLoadHints(QLibrary::ExportExternalSymbolsHint) ;
		loader.setFileName(pluginsDir.absoluteFilePath(fileName));
		
		// Convert filename for outputing		
		std::string filename(fileName.toStdString()) ;

		if(!loader.load())
		{
#ifdef DEBUG
			std::cout << "<<DEBUG>> " << loader.errorString().toStdString() << std::endl ;
#endif
			continue ;
		}

		QObject *plugin = loader.instance();
		if (plugin != NULL) 
		{
#ifdef DEBUG
			std::cout << "<<DEBUG>> loading plugin " << filename << std::endl ;
#endif
			if(dynamic_cast<function*>(plugin) != NULL)
			{
#ifdef DEBUG
				std::cout << "<<DEBUG>>  -> it is a function" << std::endl ;
#endif
				_functions.insert(std::pair<std::string, function*>(filename, dynamic_cast<function*>(plugin))) ;
			}
			
			if(dynamic_cast<data*>(plugin) != NULL)
			{
#ifdef DEBUG
				std::cout << "<<DEBUG>>  -> it is a data loader" << std::endl ;
#endif
				_datas.insert(std::pair<std::string, data*>(filename, dynamic_cast<data*>(plugin))) ;
			}
			if(dynamic_cast<fitter*>(plugin) != NULL)
			{
#ifdef DEBUG
				std::cout << "<<DEBUG>>  -> it is a fitter" << std::endl ;
#endif
				_fitters.insert(std::pair<std::string, fitter*>(filename, dynamic_cast<fitter*>(plugin))) ;
			}
		}
		else
		{
#ifdef DEBUG
			std::cout << "<<DEBUG>> " << loader.errorString().toStdString() << std::endl ;
#endif
		}
	}
	*/
}

// Get instances of the function, the data and the
// fitter
//
function* plugins_manager::get_function()
{
#ifdef USING_STATIC
    return new rational_function() ;
#else
    if(_functions.empty())
	{
		return new rational_function() ;
	}
	else
	{
		return _functions.begin()->second ;
	}
#endif

}
data* plugins_manager::get_data()
{
	//if(_datas.empty())
	{
#ifdef DEBUG
        std::cout << "<<DEBUG>>  using vertical segment data loader" << std::endl ;
#endif
		return new vertical_segment() ;
	}/*
	else
	{
#ifdef DEBUG
		std::cout << "<<DEBUG>>  using \"" << _datas.begin()->first << "\" data loader" << std::endl ;
#endif
		return _datas.begin()->second ;
	}*/
}
fitter* plugins_manager::get_fitter()
{
#ifdef USING_STATIC
    return NULL;
#else
	if(_fitters.empty())
	{
		return NULL ;
	}
	else
	{
#ifdef DEBUG
		std::cout << "<<DEBUG>>  using \"" <<  _fitters.begin()->first << "\"" << std::endl ;
#endif
		return _fitters.begin()->second ;
	}
#endif
}

arguments create_arguments(const std::string& n)
{
    std::vector<std::string> cmd_vec;
    std::stringstream stream(n);
    while(stream.good())
    {
        std::string temp;
        stream >> temp;

        cmd_vec.push_back(temp);
    }
    int argc = cmd_vec.size();
    char* argv[argc];
    for(int i=0; i<argc; ++i)
    {
        argv[i] = &cmd_vec[i][0];
    }

    arguments current_args(argc, argv);
    return current_args;
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

        //! \todo create a <em>compound</em> class to store multiple
        //! functions in it.

        //! For each args_vec element, create a function object and add
        //! it to the compound one.

        std::string n("--func ");
        n.append(args_vec[0]);
        func = get_function(create_arguments(n));

        //! return the compound class

    }
    else
    {
        std::string filename = args["func"];
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
        std::cout << "<<DEBUG>> multiplying by a Fresnel term" << std::endl;

        std::string n("--func ");
        if(args.is_vec("fresnel"))
        {
            std::string temp = args["fresnel"];
            n.append(temp.substr(1, temp.size()-2));
        }
        else
        {
            n.append(args["fresnel"]);
        }

        //fresnel* func_fres = dynamic_cast<fresnel*>(get_function(create_arguments(n)));
        //func_fres->setBase(func)
        //func = dynamic_cast<function*>(func_fres);
    }

    return func;
}
data* plugins_manager::get_data(const std::string& n)
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
fitter* plugins_manager::get_fitter(const std::string& n)
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

void plugins_manager::check_compatibility(data*& d, function*& f,
                                          const arguments& args)
{
	if(d->parametrization() == params::UNKNOWN_INPUT)
	{
		std::cout << "<<WARNING>> unknown parametrization for data" << std::endl;
	}

	if(f->parametrization() == params::UNKNOWN_INPUT)
	{
		std::cout << "<<DEBUG>> function will take the parametrization of the data" << std::endl;
		f->setParametrization(d->parametrization());
	}
	else if(d->parametrization() != f->parametrization())
	{
		std::cout << "<<INFO>> has to change the parametrization of the input data" << std::endl;
		data_params* dd = new data_params(d, f->parametrization());
		d = dd ;
	}
	else
	{
		std::cout << "<<DEBUG>> no change was made to the parametrization" << std::endl;
	}

	if(f->dimY() != d->dimY())
	{
		std::cout << "<<WARNING>> the data and the function have different Y dimensions" << std::endl;
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
#else
#include <unistd.h>
size_t plugins_manager::get_system_memory()
{
	long pages = sysconf(_SC_PHYS_PAGES);
	long page_size = sysconf(_SC_PAGE_SIZE);
	return pages * page_size;
}
#endif
