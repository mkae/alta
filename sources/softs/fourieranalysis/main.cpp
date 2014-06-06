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
		std::cout << "<<HELP>> fouriertransform --input data.file --output gnuplot.file --loader loader.so --in [input vector] --raw gnuplot.file" << std::endl ;
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

	ptr<data> d = NULL ;
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
		
		vec in(3), out(3) ;
		vec arg_in = args.get_vec("in", 2);
		double theta_in = (double)arg_in[0];
		double phi_in   = (double)arg_in[1];
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
				xdata[j] = theta;

				data[j]  = v[0] ;

				/*
				file << phi << "\t" << theta << "\t" ;
				for(int u=0; u<d->dimY(); ++u)
					file << v[u] << "\t" ;
			
				file << std::endl ;
				*/
			}

		// Create the MATLAB defintion of objects
		mxArray *d, *x, *f;
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
		engEvalString(ep, "d = d ./ normd;");
		engEvalString(ep, "ff = fftshift(fft(d));");
		engEvalString(ep, "f = sqrt(real(ff).^2 + imag(ff).^2);");
		engEvalString(ep, "normf = sum(f) / max(size(f));");
		engEvalString(ep, "f = f ./ normf;");
//		engEvalString(ep, "plot(x, f, x, d)");
//		std::cout << output << std::endl ;

		f = engGetVariable(ep, "f");
		std::ofstream out_fourier(args["output"].c_str(), std::ios_base::trunc);
		double* val = (double*)mxGetData(f) ;
		for(int i=0; i<N; ++i)
		{
			out_fourier << i - N/2 << "\t" << val[i] << std::endl ;
		}
		out_fourier.close();

		if(args.is_defined("raw"))
		{
			std::ofstream out_raw(args["raw"].c_str(), std::ios_base::trunc);
			for(int i=0; i<N; ++i)
			{
				out_raw << xdata[i] << "\t" << data[i] << std::endl ;
			}
			out_raw.close();
		}

		mxDestroyArray(x);
		mxDestroyArray(d);
		engClose(ep);
	}	
	else
	{
		std::cerr << "<<ERROR>> cannot export data" << std::endl ;
	}

	return 0 ;
}
