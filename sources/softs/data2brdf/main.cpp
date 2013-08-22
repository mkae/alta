#include <core/args.h>
#include <core/data.h>
#include <core/function.h>
#include <core/fitter.h>
#include <core/plugins_manager.h>

#include <QPluginLoader>
#include <QtPlugin>
#include <QCoreApplication>
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
	QCoreApplication app(argc, argv, false);
	arguments args(argc, argv) ;

	plugins_manager manager(args) ;
	fitter* fit = manager.get_fitter(args["fitter"]) ;
	if(fit == NULL)
	{
		fit = manager.get_fitter() ;
	}

	if(args.is_defined("available_params"))
	{
		params::print_input_params();
		return 0;
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

        function* f = plugins_manager::get_function(args);
		data*     d = plugins_manager::get_data(args["data"]);
		d->load(args["input"], args);

        if(f == NULL || d == NULL)
        {
            std::cerr << "<<ERROR>> no function or data object correctly defined" << std::endl;
            return 1;
        }

		// Check the compatibility between the data and the function
		plugins_manager::check_compatibility(d, f, args);


		// Start a timer
		QTime time ;
		time.start() ;

		// Fit the data
		bool is_fitted = fit->fit_data(d, f, args) ;

		// Get the fitting duration
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
			  double L2   = f->L2_distance(d);
			  double Linf = f->Linf_distance(d);
			  std::cout << "<<INFO>> L2   distance to data = " << L2   << std::endl;
			  std::cout << "<<INFO>> Linf distance to data = " << Linf << std::endl;

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
			file.close();
			//*/
			std::string error_filename = args["output"].substr(0,n);
			error_filename.append(".errorplot") ;
			file.open(error_filename.c_str(), std::ios_base::trunc);
			for(int i=0; i<d->size(); ++i)
			{
				vec v = d->get(i) ;
				vec y1(d->dimY()) ;
				for(int k=0; k<d->dimY(); ++k) { y1[k] = v[d->dimX() + k] ; }

				vec y2 = f->value(v) ;
				for(int u=0; u<d->dimX(); ++u)
					file << v[u] << "\t" ;

				for(int u=0; u<d->dimY(); ++u)
					file << y2[u]-y1[u] << "\t" ;

				file << std::endl ;
			}
			file.close();

			std::string linerror_filename = args["output"].substr(0,n);
			linerror_filename.append(".linearerrorplot") ;
			file.open(linerror_filename.c_str(), std::ios_base::trunc);
			for(int i=0; i<d->size(); ++i)
			{
				vec v = d->get(i) ;

				vec y1(d->dimY()) ;
				for(int k=0; k<d->dimY(); ++k) { y1[k] = 0.5*(v[d->dimX() + k] +v[d->dimX()+d->dimY() + k]); }
				vec y2 = f->value(v) ;

				file << i << "\t" ;
				for(int u=0; u<d->dimY(); ++u)
					file << y2[u]-y1[u] << "\t" ;

				file << std::endl ;
			}
			file.close();
			//*/
#endif
			return 0;
		}
		else
		{
			std::cout << "<<ERROR>> unable to fit the data" << std::endl ;
			return 1;
		}

	}	
	/*	else
		{
		std::cout << "<<ERROR>> not enough plugin defined" << std::endl ;
		}
		*/

	return 0 ;
}
