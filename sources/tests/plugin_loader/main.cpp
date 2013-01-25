#include <core/args.h>
#include <core/data.h>
#include <core/function.h>
#include <core/fitter.h>

#include <QPluginLoader>
#include <QtPlugin>
#include <QApplication>
#include <QDir>

#include <iostream>
#include <vector>
#include <iostream>
#include <fstream>
#include <limits>

int main(int argc, char** argv)
{
	QApplication app(argc, argv);


	std::vector<function*> functions ;
	std::vector<data*>     datas ;
	std::vector<fitter*>   fitters ;

	QDir pluginsDir = QDir(app.applicationDirPath());
	foreach (QString fileName, pluginsDir.entryList(QDir::Files)) 
	{
		QPluginLoader loader(pluginsDir.absoluteFilePath(fileName));

		QObject *plugin = loader.instance();
		if (plugin != nullptr) 
		{
#ifdef DEBUG
			std::cout << "<<DEBUG>> loading plugin " << fileName.toStdString() << std::endl ;
#endif
			if(dynamic_cast<function*>(plugin) != nullptr)
			{
#ifdef DEBUG
				std::cout << "<<DEBUG>>  -> it is a function" << std::endl ;
#endif
				functions.push_back(dynamic_cast<function*>(plugin)) ;
			}
			else if(dynamic_cast<data*>(plugin) != nullptr)
			{
#ifdef DEBUG
				std::cout << "<<DEBUG>>  -> it is a data loader" << std::endl ;
#endif
				datas.push_back(dynamic_cast<data*>(plugin)) ;
			}
			else if(dynamic_cast<fitter*>(plugin) != nullptr)
			{
#ifdef DEBUG
				std::cout << "<<DEBUG>>  -> it is a fitter" << std::endl ;
#endif
				fitters.push_back(dynamic_cast<fitter*>(plugin)) ;
			}

		}
	}

	arguments args(argc, argv) ;
	if(! args.is_defined("input")) {
		std::cerr << "<<ERROR>> the input filename is not defined" << std::endl ;
		return 1 ;
	}
	if(! args.is_defined("output")) {
		std::cerr << "<<ERROR>> the output filename is not defined" << std::endl ;
		return 1 ;
	}

	if(fitters.size() > 0 && datas.size() > 0 && functions.size() > 0)
	{
		fitters[0]->set_parameters(args) ;

		function* f = functions[0] ;
		data*     d = datas[0] ;
		if(args.is_defined("min") && args.is_defined("max"))
		{
			d->load(args["input"], args.get_float("min", 0.0f), args.get_float("max", 1.0f));
		}
		else if(args.is_defined("min") && !args.is_defined("max"))
		{
			d->load(args["input"], args.get_float("min", 0.0f), std::numeric_limits<double>::max());
		}
		else if(args.is_defined("min") && !args.is_defined("max"))
		{
			d->load(args["input"], -std::numeric_limits<double>::max(), args.get_float("min", 0.0f));
		}
		else
		{
			d->load(args["input"]);
		}

		bool is_fitted = fitters[0]->fit_data(d, f) ;

		// Display the result
		if(is_fitted)
		{
//			std::cout << *f << std::endl ;

			//*
			std::ofstream file(args["output"], std::ios_base::trunc);
			const double dt = (d->max() - d->min()) / 100.0f ;
			for(double x=d->min(); x<=d->max(); x+=dt)
			{
				file << x << "\t" << (*f)(x) << std::endl ;
			}
			//*/
		}
		else
		{
			std::cout << "<<ERROR>> unable to fit the data" << std::endl ;
		}

	}	
	else
	{
		std::cout << "<<ERROR>> not enough plugin defined" << std::endl ;
	}


	return 0 ;
}
