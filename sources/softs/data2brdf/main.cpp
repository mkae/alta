#include <core/args.h>
#include <core/data.h>
#include <core/function.h>
#include <core/fitter.h>
#include <core/plugins_manager.h>

#include <iostream>
#include <vector>
#include <iostream>
#include <fstream>
#include <limits>
#include <cstdlib>

int main(int argc, char** argv)
{
    arguments args(argc, argv) ;

    fitter* fit = plugins_manager::get_fitter(args["fitter"]) ;
    if(fit == NULL)
    {
        fit = plugins_manager::get_fitter() ;
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
    if(fit != NULL)
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

		  if(d->size() == 0)
		  {
			  std::cerr << "<<ERROR>> no data loaded, please check you input file" << std::endl;
			  return 1;
		  }

        // Check the compatibility between the data and the function
        plugins_manager::check_compatibility(d, f, args);


        // Start a timer
        timer time ;
        time.start() ;

        // Fit the data
        bool is_fitted = fit->fit_data(d, f, args) ;

        // Get the fitting duration
        time.stop();

        // Display the result
        if(is_fitted)
        {
            std::cout << "<<INFO>> total time: " << time << std::endl ;

            double L2   = f->L2_distance(d);
            double Linf = f->Linf_distance(d);
            std::cout << "<<INFO>> L2   distance to data = " << L2   << std::endl;
            std::cout << "<<INFO>> Linf distance to data = " << Linf << std::endl;

            // Export the L2 and Linf values to the command line
            std::stringstream L2string, Linfstring;
            L2string << L2; Linfstring << Linf;
            args.update("L2",   L2string.str());
            args.update("Linf", Linfstring.str());

            f->save(args["output"], args) ;
            return 0;
        }
        else
        {
            std::cout << "<<ERROR>> unable to fit the data" << std::endl ;
            return 1;
        }

    }
    else
    {
        std::cout << "<<ERROR>> no fitter loaded, please check your command line arguments" << std::endl ;
    }

    return 0 ;
}
