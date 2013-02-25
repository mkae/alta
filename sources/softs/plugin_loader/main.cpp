#include <core/args.h>
#include <core/data.h>
#include <core/function.h>
#include <core/fitter.h>
#include <core/plugins_manager.h>

#include <QPluginLoader>
#include <QtPlugin>
#include <QApplication>
#include <QDir>
#include <QTime>

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

	plugins_manager manager(args) ;
	fitter* fit = manager.get_fitter(args["fitter"]) ;
	if(fit == NULL)
	{
		fit = manager.get_fitter() ;
	}

	if(! args.is_defined("input")) {
		std::cerr << "<<ERROR>> the input filename is not defined" << std::endl ;
		return 1 ;
	}
	if(! args.is_defined("output")) {
		std::cerr << "<<ERROR>> the output filename is not defined" << std::endl ;
		return 1 ;
	}

//	if(fitters.size() > 0 && datas.size() > 0 && functions.size() > 0)
	{
		fit->set_parameters(args) ;
		function* f = NULL;
		if(args.is_defined("func"))
		{
			std::cout << "<<INFO>> Using plugin function \"" << args["func"] << "\"" << std::endl ;
			f = manager.get_function(args["func"]) ;
		}
		else
		{
			f = fit->provide_function() ;
		}
		data*     d = fit->provide_data() ;
		d->load(args["input"], args);

		QTime time ;
		time.start() ;
		bool is_fitted = fit->fit_data(d, f) ;
		int msec = time.elapsed() ;
		int sec  = (msec / 1000) % 60 ;
		int min  = (msec / 60000) % 60 ;
		int hour = (msec / 3600000) ;
		

		// Display the result
		if(is_fitted)
		{
			std::cout << "<<INFO>> total time: " << hour << "h " << min << "m " << sec << "s" << std::endl ;
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
#ifdef OLD // use brdf2gnuplot
			size_t n = args["output"].find('.') ;
			std::string gnuplot_filename = args["output"].substr(0,n); 
			gnuplot_filename.append(".gnuplot") ;
/*	
			f->save_gnuplot(gnuplot_filename, d, args);				
/*/
			std::ofstream file(gnuplot_filename.c_str(), std::ios_base::trunc);
			for(int i=0; i<d->size(); ++i)
			{
				vec v = d->get(i) ;
//				vec y1(d->dimY()) ;
//				for(int k=0; k<d->dimY(); ++k) { y1[k] = v[d->dimX() + k] ; }

				vec y2 = f->value(v) ;
				for(int u=0; u<d->dimX(); ++u)
					file << v[u] << "\t" ;
					
				for(int u=0; u<d->dimY(); ++u)
					file << y2[u] << "\t" ;
					
				file << std::endl ;
			}	
//*/
			std::string error_filename = args["output"].substr(0,n); 
			error_filename.append(".errorplot") ;
			std::ofstream efile(error_filename.c_str(), std::ios_base::trunc);
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
#endif
		}
		else
		{
			std::cout << "<<ERROR>> unable to fit the data" << std::endl ;
		}

	}	
/*	else
	{
		std::cout << "<<ERROR>> not enough plugin defined" << std::endl ;
	}
*/

	return 0 ;
}
