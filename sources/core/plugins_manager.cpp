#include "plugins_manager.h"
#include "rational_function.h"
#include "vertical_segment.h"

#include <QCoreApplication>
#include <QPluginLoader>
#include <QtPlugin>
#include <QLibrary>
#include <QDir>

// Create the object, parse the argument and load
// all the plugins
//
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
}

#define USING_STATIC
#ifdef USING_STATIC
typedef function* (*FunctionPrototype)();
typedef fitter*   (*FitterPrototype)();
typedef data*     (*DataPrototype)();
#endif

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
data* plugins_manager::get_data()     const 
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
fitter* plugins_manager::get_fitter()   const 
{
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
}


// Get instances of the function, the data and the
// fitter, select one based on the name. Return null
// if no one exist.
//
function* plugins_manager::get_function(const std::string& n)
{
#ifdef USING_STATIC
    FunctionPrototype myFunction = (FunctionPrototype) QLibrary::resolve(QString(n.c_str()), "_Z16provide_functionv");

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
data* plugins_manager::get_data(const std::string& n)     const 
{
	if(_datas.count(n) == 0)
	{
#ifdef DEBUG
		std::cout << "<<DEBUG>>  using vertical segment data loader" << std::endl ;
#endif
		return new vertical_segment() ;
	}
	else
	{
#ifdef DEBUG
		std::cout << "<<DEBUG>>  using \"" << n << "\" data loader" << std::endl ;
#endif
		return _datas.find(n)->second ;
	}
}
fitter* plugins_manager::get_fitter(const std::string& n)   const 
{
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
