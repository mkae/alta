#include <core/args.h>
#include <core/function.h>
#include <core/plugins_manager.h>

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
	QCoreApplication app(argc, argv);
	arguments args(argc, argv) ;

	if(! args.is_defined("input")) {
		std::cerr << "<<ERROR>> the input filename is not defined" << std::endl ;
		return 1 ;
	}
	if(! args.is_defined("output")) {
		std::cerr << "<<ERROR>> the output filename is not defined" << std::endl ;
		return 1 ;
	}

    // Load the function
    function* f = plugins_manager::get_function(args["input"]);

    // Save it
    f->save(args["output"], args) ;

	return 0 ;
}
