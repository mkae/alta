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
		if (plugin != nullptr) {

			std::cout << "Loading plugin " << fileName.toStdString() << std::endl ;
			if(dynamic_cast<function*>(plugin) != nullptr)
			{
				std::cout << "  -> it is a function" << std::endl ;
				functions.push_back(dynamic_cast<function*>(plugin)) ;
			}
			else if(dynamic_cast<data*>(plugin) != nullptr)
			{
				std::cout << "  -> it is a data loader" << std::endl ;
				datas.push_back(dynamic_cast<data*>(plugin)) ;
			}
			else if(dynamic_cast<fitter*>(plugin) != nullptr)
			{
				std::cout << "  -> it is a fitter" << std::endl ;
				fitters.push_back(dynamic_cast<fitter*>(plugin)) ;
			}

		}
	}

	arguments args(argc, argv) ;

	if(fitters.size() > 0 && datas.size() > 0 && functions.size() > 0)
	{
		fitters[0]->set_parameters(args) ;

		function* f = functions[0] ;
		data*     d = datas[0] ;
		bool is_fitted = fitters[0]->fit_data(d, f) ;

		// Display the result
		if(is_fitted)
		{
			std::cout << f << std::endl ;

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
