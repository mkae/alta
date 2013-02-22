#include <core/args.h>
#include <core/data.h>
#include <core/function.h>
#include <core/fitter.h>
#include <core/plugins_manager.h>

#include <QApplication>

#include <iostream>
#include <vector>
#include <iostream>
#include <fstream>
#include <limits>
#include <cstdlib>
#include <cmath>
#include <engine.h>

#define BUFFER_SIZE 10000

int main(int argc, char** argv)
{
	QApplication app(argc, argv, false);
	arguments args(argc, argv) ;

	plugins_manager manager(args) ;

	if(args.is_defined("help")) {
		std::cout << "<<HELP>> data2gnuplot --input data.file --output gnuplot.file --loader loader.so" << std::endl ;
		std::cout << " - input and output are mandatory parameters" << std::endl ;
		return 0 ;
	}

	if(! args.is_defined("input")) {
		std::cerr << "<<ERROR>> the input filename is not defined" << std::endl ;
		return 1 ;
	}
	if(! args.is_defined("output")) {
		std::cerr << "<<ERROR>> the output filename is not defined" << std::endl ;
		return 1 ;
	}

	data* d = NULL ;
	d = manager.get_data(args["loader"]) ;
	d->load(args["input"]);

	// Create output file
	std::ofstream file(args["output"].c_str(), std::ios_base::trunc);

	if(d != NULL)
	{
		// Create matlab engine
		Engine *ep;
		if (!(ep = engOpen(""))) 
		{
			std::cerr << "<ERROR>> can't start MATLAB engine" << std::endl ;
			return false ;
		}


		std::cout << "<<INFO>> will export " << d->size() << " elements" << std::endl ;
	
		double theta_in = (double)args.get_float("theta", 0.0f);
		double phi_in   = (double)args.get_float("phi", 0.0f);
		vec in(3), out(3) ;
		in[0] = cos(phi_in)*sin(theta_in);
		in[1] = sin(phi_in)*sin(theta_in);
		in[2] = cos(theta_in);

		const int N = 100;
		double* data  = new double[N];
		double* xdata = new double[N];
//		for(int i=0; i<N; ++i)
		int i=0;
		for(int j=0; j<N; ++j)
			{
				double phi   = (i-N/2) * M_PI / (N-1) * 2;
				double theta = (j-N/2) * M_PI / (N-1) ;
				out[0] = cos(phi)*sin(theta);
				out[1] = sin(phi)*sin(theta);
				out[2] = cos(theta);
				vec v = d->value(in, out) ;
				xdata[j] = j-N/2;
				data[j]  = v[0] ;

				/*
				file << phi << "\t" << theta << "\t" ;
				for(int u=0; u<d->dimY(); ++u)
					file << v[u] << "\t" ;
			
				file << std::endl ;
				*/
			}

		// Create the MATLAB defintion of objects
		// MATLAB defines a quad prog as
		//   1/2 x' H x + f' x with A x <= b 
		//
		mxArray *d, *x;
		d = mxCreateDoubleMatrix(  N, 1, mxREAL);
		x = mxCreateDoubleMatrix(  N, 1, mxREAL);
		memcpy((void *)mxGetPr(d), (void *) data,  N*sizeof(double));
		memcpy((void *)mxGetPr(x), (void *) xdata,  N*sizeof(double));

		engPutVariable(ep, "x", x);
		engPutVariable(ep, "d", d);

		char* output = new char[BUFFER_SIZE+1];
		output[BUFFER_SIZE] = '\0';
		engOutputBuffer(ep, output, BUFFER_SIZE) ;
		engEvalString(ep, "normd = sum(d) / max(size(d));");
		std::cout << output << std::endl ;
		engEvalString(ep, "d = d ./ normd;");
		std::cout << output << std::endl ;
		engEvalString(ep, "f = real(fftshift(fft(d)));");
		std::cout << output << std::endl ;
		engEvalString(ep, "normf = sum(f) / max(size(f));");
		std::cout << output << std::endl ;
		engEvalString(ep, "f = f ./ normf;");
		std::cout << output << std::endl ;
		engEvalString(ep, "plot(x, f, x, d)");
		std::cout << output << std::endl ;

		char c;
		std::cin >> c ;

		engClose(ep);
	}	
	else
	{
		std::cerr << "<<ERROR>> cannot export data" << std::endl ;
	}

	return 0 ;
}
