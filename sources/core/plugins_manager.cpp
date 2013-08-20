#include "plugins_manager.h"
#include "rational_function.h"
#include "vertical_segment.h"

#include <QCoreApplication>
/*
#include <QPluginLoader>
#include <QtPlugin>
#include <QLibrary>
*/
#include <QDir>


#ifdef WIN32
#else
    #include <stdio.h>
    #include <dlfcn.h>
#endif

//! \brief Open a dynamic library file (.so or .dll) and extract the associated
//! provide function. The template argument is used to cast the library to a
//! specific type.
template<typename T> T open_library(const std::string& filename, const char* function)
{
#ifdef WIN32
    return NULL;
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
        return NULL;
    }
#endif
}


// Create the object, parse the argument and load all the plugins
plugins_manager::plugins_manager(const arguments& args) 
{

	QDir pluginsDir;
	if(args.is_defined("plugins"))
	{
		pluginsDir = QDir(args["plugins"].c_str()) ;
	}
	else if(QCoreApplication::instance() != NULL)
	{
		pluginsDir = QDir(QCoreApplication::instance()->applicationDirPath());
	}

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

// Get instances of the function, the data and the
// fitter, select one based on the name. Return null
// if no one exist.
//
function* plugins_manager::get_function(const std::string& n)
{
    if(n.empty())
    {
#ifdef DEBUG
        std::cout << "<<DEBUG>> no function plugin specified, returning a rational function" << std::endl;
#endif
        return new rational_function();
    }

#ifdef USING_STATIC
	 std::string file;
	 if(n[0] == '.')
	 {
		 file = n.substr(1, n.size()-1);
	 }
	 else
	 {
		 file = n;
	 }

    QString path = QDir::currentPath() + QString(file.c_str()) ;
	 
	 FunctionPrototype myFunction = open_library<FunctionPrototype>(path.toStdString(), "provide_function");
    if(myFunction != NULL)
    {
#ifdef DEBUG
        std::cout << "<<DEBUG>> using function provider in file \"" << n << "\"" << std::endl;
#endif
        return myFunction();
    }
    else
    {
        std::cerr << "<<ERROR>> no function provider found in file \"" << n << "\"" << std::endl;
        return new rational_function() ;
    }
#else

	if(_functions.count(n) == 0)
	{
		return new rational_function() ;
	}
	else
	{
		return _functions.find(n)->second ;
	}
#endif
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

#ifdef USING_STATIC
	 std::string file;
	 if(n[0] == '.')
	 {
		 file = n.substr(1, n.size()-1);
	 }
	 else
	 {
		 file = n;
	 }

    QString path = QDir::currentPath() + QString(file.c_str()) ;
	 
	 DataPrototype myData = open_library<DataPrototype>(path.toStdString(), "provide_data");
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
#else

    if(_functions.count(n) == 0)
    {
        return new vertical_segment() ;
    }
    else
    {
        return _datas.find(n)->second ;
    }
#endif
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

#ifdef USING_STATIC
	 std::string file;
	 if(n[0] == '.')
	 {
		 file = n.substr(1, n.size()-1);
	 }
	 else
	 {
		 file = n;
	 }

    QString path = QDir::currentPath() + QString(file.c_str()) ;
    
	 FitterPrototype myFitter = open_library<FitterPrototype>(path.toStdString(), "provide_fitter");

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
#else
	if(_fitters.count(n) == 0)
	{
		return NULL ;
	}
	else
	{
#ifdef DEBUG
		std::cout << "<<DEBUG>>  using \"" <<  n << "\"" << std::endl ;
#endif
		return _fitters.find(n)->second ;
	}
#endif
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
