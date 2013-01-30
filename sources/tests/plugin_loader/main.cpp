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
#include <cstdlib>

int main(int argc, char** argv)
{
//	QCoreApplication::addLibraryPath() ;
//	QCoreApplication::addLibraryPath(QString("/home/belcour/Projects/alta/sources/tests/plugin_loader/")) ;

	QApplication app(argc, argv, false);
	arguments args(argc, argv) ;

	std::vector<function*> functions ;
	std::vector<data*>     datas ;
	std::vector<fitter*>   fitters ;

	QDir pluginsDir = QDir(app.applicationDirPath());
	if(args.is_defined("plugins"))
	{
		pluginsDir = QDir(args["plugins"].c_str()) ;
	}

	std::vector<QString> plugins ;

	foreach (QString fileName, pluginsDir.entryList(QDir::Files)) 
	{
		QPluginLoader loader ;
		loader.setLoadHints(QLibrary::ExportExternalSymbolsHint) ;
		loader.setFileName(pluginsDir.absoluteFilePath(fileName));


		QObject *plugin = loader.instance();
		if (plugin != NULL) 
		{
#ifdef DEBUG
			std::cout << "<<DEBUG>> loading plugin " << fileName.toStdString() << std::endl ;
#endif
			if(dynamic_cast<function*>(plugin) != NULL)
			{
#ifdef DEBUG
				std::cout << "<<DEBUG>>  -> it is a function" << std::endl ;
#endif
				functions.push_back(dynamic_cast<function*>(plugin)) ;
			}
			
			if(dynamic_cast<data*>(plugin) != NULL)
			{
#ifdef DEBUG
				std::cout << "<<DEBUG>>  -> it is a data loader" << std::endl ;
#endif
				datas.push_back(dynamic_cast<data*>(plugin)) ;
			}
			
			if(dynamic_cast<fitter*>(plugin) != NULL)
			{
#ifdef DEBUG
				std::cout << "<<DEBUG>>  -> it is a fitter" << std::endl ;
#endif
				fitters.push_back(dynamic_cast<fitter*>(plugin)) ;
			}

		}
		else
		{
#ifdef DEBUG
			std::cout << "<<DEBUG>> " << loader.errorString().toStdString() << std::endl ;
#endif
			plugins.push_back(fileName) ;
		}
	}


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
		d->load(args["input"], args);

		bool is_fitted = fitters[0]->fit_data(d, f) ;

		// Display the result
		if(is_fitted)
		{
/*
			vec min, max ;
			min.assign(2, args.get_float("min", 0.0f)) ;
			max.assign(2, args.get_float("max", 1.5f)) ;

			int nb_samples = args.get_int("nb_samples", 100) ;
			double dt = (max[0]-min[0]) / nb_samples ;

			std::ofstream file(args["output"].c_str(), std::ios_base::trunc);
			for(double x=min[0]; x<=max[0]; x+=dt)
			{
				vec vx ; for(int i=0;i<2; ++i) { vx.push_back(x) ; }
				file << x << "\t" << f->value(vx)[0] << std::endl ;
				std::cout << x << "\t" << f->value(vx)[0] << std::endl ;
			}
/*/
			f->save(args["output"], args) ;

			std::ofstream file("temp.gnuplot", std::ios_base::trunc);
			for(int i=0; i<d->size(); ++i)
			{
				vec v = d->get(i) ;
				vec y1(d->dimY()) ;
				for(int k=0; k<d->dimY(); ++k) { y1[k] = v[d->dimX() + k] ; }

				vec y2 = f->value(v) ;
				for(int u=0; u<d->dimX(); ++u)
					file << v[u] << "\t" ;
					
				for(int u=0; u<d->dimY(); ++u)
					file << y2[u] << "\t" ;
					
				file << std::endl ;
			}	

			std::ofstream efile("error.gnuplot", std::ios_base::trunc);
			for(int i=0; i<d->size(); ++i)
			{
				vec v = d->get(i) ;
				vec y1(d->dimY()) ;
				for(int k=0; k<d->dimY(); ++k) { y1[k] = v[d->dimX() + k] ; }

				vec y2 = f->value(v) ;
				for(int u=0; u<d->dimX(); ++u)
					efile << v[u] << "\t" ;
					
				for(int u=0; u<d->dimY(); ++u)
					efile << y2[u]-y1[u] << "\t" ;
					
				efile << std::endl ;
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
