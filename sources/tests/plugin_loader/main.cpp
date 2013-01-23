#include <core/args.h>
#include <core/data.h>
#include <core/function.h>
#include <core/fitter.h>

#include <QPluginLoader>
#include <QtPlugin>
#include <QApplication>
#include <QDir>

#include <iostream>

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
		if (plugin) {

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


	std::cout << "Done" << std::endl ;

	return 0 ;
}
