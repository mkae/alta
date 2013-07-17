#include "rational_fitter.h"

#include <Eigen/Dense>
#include <Eigen/SVD>

#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>

#include <QTime>

#define BUFFER_SIZE 10000

ALTA_DLL_EXPORT fitter* provide_fitter()
{
    return new rational_fitter_dca();
}

data* rational_fitter_dca::provide_data() const
{
	return new vertical_segment() ;
}

function* rational_fitter_dca::provide_function() const 
{
	return new rational_function() ;
}

rational_fitter_dca::rational_fitter_dca() : QObject()
{
}
rational_fitter_dca::~rational_fitter_dca() 
{
}

bool rational_fitter_dca::fit_data(const data* dat, function* fit, const arguments &args)
{
	rational_function* r = dynamic_cast<rational_function*>(fit) ;
	const vertical_segment* d = dynamic_cast<const vertical_segment*>(dat) ;
	if(r == NULL || d == NULL)
	{
		std::cerr << "<<ERROR>> not passing the correct class to the fitter" << std::endl ;
		return false ;
	}

	// Create matlab engine
	if (!(ep = engOpen(""))) 
	{
		std::cerr << "<ERROR>> can't start MATLAB engine" << std::endl ;
		return false ;
	}


	// I need to set the dimension of the resulting function to be equal
	// to the dimension of my fitting problem
	r->setDimX(d->dimX()) ;
	r->setDimY(d->dimY()) ;
	r->setMin(d->min()) ;
	r->setMax(d->max()) ;


	QTime time ;
	time.start() ;

	if(fit_data(d, _max_np, _max_nq, r, args))
	{
		int msec = time.elapsed() ;
		int sec  = (msec / 1000) % 60 ;
		int min  = (msec / 60000) % 60 ;
		int hour = (msec / 3600000) ;
		std::cout << "<<INFO>> got a fit using np = " << _max_np << " & nq =  " << _max_nq << "      " << std::endl ;
		std::cout << "<<INFO>> it took " << hour << "h " << min << "m " << sec << "s" << std::endl ;

		return true ;
	}

	engClose(ep); 
	return false ;
}

void rational_fitter_dca::set_parameters(const arguments& args)
{
	_max_np = args.get_float("np", 10) ;
	_max_nq = args.get_float("nq", 10) ;
}

double distance(const rational_function* f, const data* d)
{
	double distance = 0.0;
	for(int i=0; i<d->size(); ++i)
	{
		vec xi = d->get(i) ;
		vec y  = f->value(xi);
		double current_d = 0.0;

		for(int j=0; j<d->dimY(); ++j)
		{
			double diff = y[j] - xi[d->dimX()+j];
			current_d += diff;
		}

		current_d = std::abs(current_d);
		if(current_d > distance)
			distance = current_d;
	}
	return distance;
}

// Bootstrap the DCA algorithm with an already done fit
void rational_fitter_dca::bootstrap(const data* d, int np, int nq, rational_function* fit, double& delta, const arguments& args)
{
	vec p(np*d->dimY());
	vec q(nq*d->dimY());
	
	if(args.is_defined("bootstrap"))
	{
		fit->load(args["bootstrap"]);
	}
	else
	{
#ifdef DEBUG
		std::cout << "<<DEBUG>> Using the constant function equals to 0 as input: not optimal" << std::endl;
#endif
		q[0] = 1.0;
		fit->update(p, q);
	}

	delta = distance(fit, d);
}

// dat is the data object, it contains all the points to fit
// np and nq are the degree of the RP to fit to the data
// y is the dimension to fit on the y-data (e.g. R, G or B for RGB signals)
// the function return a ration BRDF function and a boolean
bool rational_fitter_dca::fit_data(const data* d, int np, int nq, rational_function* r, const arguments& args)
{
	// Size of the problem
	int N  = np+nq+1 ;
	int M  = d->size() ;
	int nY = d->dimY();	

	// Bootstrap the delta and rational function using the Papamarkos
	// algorithm.
	double delta = 0.0;
	bootstrap(d, np, nq, r, delta, args);
#ifdef DEBUG
	std::cout << "<<DEBUG>> delta value after boostrap: " << delta << std::endl;
	std::cout << "<<DEBUG>> r: " << *r << std::endl;
#endif

	// Create the MATLAB defintion of objects
	// MATLAB defines a linear prog as
	//   min f' x with A x <= b 
	//
	mxArray *f, *A, *b, *x, *flag, *u, *l;
	f = mxCreateDoubleMatrix( N*nY,      1, mxREAL);
	A = mxCreateDoubleMatrix( 4*M*nY, N*nY, mxREAL);
	b = mxCreateDoubleMatrix( 4*M*nY,    1, mxREAL);
	u = mxCreateDoubleMatrix( N*nY,      1, mxREAL);
	l = mxCreateDoubleMatrix( N*nY,      1, mxREAL);

	// Matrices of the problem in Eigen format
	Eigen::VectorXd g (nY*N) ;
	Eigen::MatrixXd CI(4*M*nY, N*nY) ;
	Eigen::VectorXd ci(4*M*nY) ;
	Eigen::VectorXd li(N*nY) ;
	Eigen::VectorXd ui(N*nY) ;

	// Updating the upper and lower values for the solution
	for(int i=0; i<N*nY; ++i)
	{
		if(i<np*nY)
		{
			ui(i) =  1.0E30; // std::numeric_limits<double>::max();
			li(i) = -1.0E30; //-std::numeric_limits<double>::max();
		}
		else if(i<(np+nq)*nY)
		{
			ui(i) =  1.0;
			li(i) = -1.0;
		}
		else
		{
			ui(i) =  1.0E30; // std::numeric_limits<double>::max();
			li(i) = -1.0E30; //-std::numeric_limits<double>::max();
		}
	}
	memcpy((void *)mxGetPr(u), (void *)ui.data(), N*nY*sizeof(double));
	memcpy((void *)mxGetPr(l), (void *)li.data(), N*nY*sizeof(double));
	engPutVariable(ep, "u", u);
	engPutVariable(ep, "l", l);

	double delta_k;
	unsigned int nb_passes = 1;

	// Loop until you get a converge solution \delta > \delta_k
	do
	{
		// delta_{k+1} = delta_{k}
		//delta_k = distance(r, d);
		delta_k = delta;

		// The function to minimize is \delta which is the last element of
		// the result vector
		for(int j=0; j<nY*(N-1); ++j)
		{
			g(j) = 0.0;
		}
		for(int y=0; y<nY; ++y)
		{
			g(nY*(N-1)+y) = 1.0;
		}

		// Each constraint (fitting interval or point add another dimension to the constraint
		// matrix
		for(int i=0; i<M; ++i)	
		{		
			//
			vec xi = d->get(i) ;

			// For each input data i \in M, the following constraints have to
			// be fulfilled:
			//   [ f_i + \delta_k] * q_i - p_i + qk_i \delta >= 0
			//   [-f_i + \delta_k] * q_i + p_i + qk_i \delta >= 0
			//                                             p >= 0
			//                                             q >  0
			//
			// The constraints matrix has the following form
			// [-p_j(x_i) ..-p_j(x_i), [ f_i+\delta_k]*q_j(x_i) ..[ fi+\delta_k]*q_j(x_i), qk(x_i)]
			// [ p_j(x_i) .. p_j(x_i), [-f_i+\delta_k]*q_j(x_i) ..[-fi+\delta_k]*q_j(x_i), qk(x_i)]
			//
			for(int j=0; j<N; ++j)
			{
				// Filling the p part
				if(j<np)
				{
					const double pi = r->p(xi, j) ;

					// Updating Eigen matrix
					for(int y=0; y<nY; ++y)
					{
						CI(2*(nY*i + y)+0, nY*j + y) =  pi ;
						CI(2*(nY*i + y)+1, nY*j + y) = -pi ;
						CI(2*M*nY + 2*(nY*i+y) + 0, nY*j+y) = -pi ;
						CI(2*M*nY + 2*(nY*i+y) + 1, nY*j+y) = 0.0 ;
					}
				}
				// Filling the q part
				else if(j<np+nq)
				{
					const double qi = r->q(xi, j-np) ;

					// Updating Eigen matrix
					for(int y=0; y<nY; ++y)
					{
						CI(2*(nY*i + y)+0, nY*j + y) = -(delta_k+xi[d->dimX()+y]) * qi ;
						CI(2*(nY*i + y)+1, nY*j + y) = -(delta_k-xi[d->dimX()+y]) * qi ;
						CI(2*M*nY + 2*(nY*i+y) + 0, nY*j + y) = 0.0 ;
						CI(2*M*nY + 2*(nY*i+y) + 1, nY*j + y) = - qi ;
					}
				}
				else
				{
					// Last column of the constraint matrix
					vec qk = r->q(xi) ;
					for(int y=0; y<nY; ++y)
					{
						CI(2*(nY*i + y)+0, nY*j + y) = -qk[y] ;
						CI(2*(nY*i + y)+1, nY*j + y) = -qk[y] ;
						CI(2*M*nY + 2*(nY*i+y) + 0, nY*j + y) = 0.0 ;
						CI(2*M*nY + 2*(nY*i+y) + 1, nY*j + y) = 0.0 ;
					}
				}
			}


			// Set the c vector
			for(int y=0; y<nY; ++y)
			{
				ci(2*(nY*i+y)+0)  = 0.0 ;
				ci(2*(nY*i+y)+1)  = 0.0 ;
				ci(2*M*nY + 2*(nY*i) + 0) = 0.0 ;
				ci(2*M*nY + 2*(nY*i) + 1) = 0.0 ;
			}
		}

		// Copy the data to matlab and execute the linear program
		//
		memcpy((void *)mxGetPr(f), (void *) g.data(), N*nY*sizeof(double));
		memcpy((void *)mxGetPr(A), (void *)CI.data(), 4*M*nY*N*nY*sizeof(double));
		memcpy((void *)mxGetPr(b), (void *)ci.data(), 4*M*nY*sizeof(double));

		engPutVariable(ep, "f", f);
		engPutVariable(ep, "A", A);
		engPutVariable(ep, "b", b);

		char* output = new char[BUFFER_SIZE+1];
		output[BUFFER_SIZE] = '\0';
		engOutputBuffer(ep, output, BUFFER_SIZE) ;

		engEvalString(ep, "[x, fval, flag] = linprog(f,A,b,[],[], l, u);");
#ifdef DEBUG
		std::cout << output << std::endl ;
#endif

		x    = engGetVariable(ep, "x") ;
		flag = engGetVariable(ep, "flag") ;

		// Update the rational function
		double* val = (double*)mxGetData(x) ;
		std::vector<double> a, b;
		for(int i=0; i<(np+nq)*nY; ++i)
		{
			if(i < np*nY)
			{
				a.push_back(val[i]) ;
			}
			else
			{
				b.push_back(val[i]) ;
			}
		}
		
		std::vector<double> tempP = r->getP();
		std::vector<double> tempQ = r->getQ();

		r->update(a, b) ;
#ifdef DEBUG
		std::cout << "<<DEBUG>> current rational function: " <<  *r << std::endl ;
#endif

		// Compute the new delta_k, the distance to the data points
		delta = distance(r, d);
		//delta = val[(np+nq)*nY];
		

		// Stopping condition if the optimization did not manage to improve the
		// L_inf norm quit !
		if(delta > delta_k)
		{

			r->update(tempP, tempQ);
			break;
		}
	
		++nb_passes;
	}while(true);

	mxDestroyArray(f);
	mxDestroyArray(A);
	mxDestroyArray(b);
	mxDestroyArray(u);
	mxDestroyArray(l);

	if(nb_passes == 1)
	{
		std::cout << "<<ERROR>> Could no optimize with respect to Linf" << std::endl;
		return false;
	}
	else
	{
		std::cout << "<<INFO>> Used " << nb_passes << " passes to optimize the solution" << std::endl;
		return true;
	}
}

Q_EXPORT_PLUGIN2(rational_fitter_dca, rational_fitter_dca)
