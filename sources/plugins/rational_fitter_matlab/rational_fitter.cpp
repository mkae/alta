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

data* rational_fitter_matlab::provide_data() const
{
	return new rational_data() ;
}

function* rational_fitter_matlab::provide_function() const 
{
	return new rational_function() ;
}

rational_fitter_matlab::rational_fitter_matlab() : QObject()
{
}
rational_fitter_matlab::~rational_fitter_matlab() 
{
}

bool rational_fitter_matlab::fit_data(const data* dat, function* fit)
{
	rational_function* r = dynamic_cast<rational_function*>(fit) ;
	const rational_data* d = dynamic_cast<const rational_data*>(dat) ;
	if(r == NULL || d == NULL || ep == NULL)
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

	std::cout << "<<INFO>> np in  [" << _min_np << ", " << _max_np 
	          << "] & nq in [" << _min_nq << ", " << _max_nq << "]" << std::endl ;

	int temp_np = _min_np, temp_nq = _min_nq ;
	while(temp_np <= _max_np || temp_nq <= _max_nq)
	{
		QTime time ;
		time.start() ;
		
		if(fit_data(d, temp_np, temp_nq, r))
		{
			int msec = time.elapsed() ;
			int sec  = (msec / 1000) % 60 ;
			int min  = (msec / 60000) % 60 ;
			int hour = (msec / 3600000) ;
			std::cout << "<<INFO>> got a fit using np = " << temp_np << " & nq =  " << temp_nq << "      " << std::endl ;
			std::cout << "<<INFO>> it took " << hour << "h " << min << "m " << sec << "s" << std::endl ;

			return true ;
		}

		std::cout << "<<INFO>> fit using np = " << temp_np << " & nq =  " << temp_nq << " failed\r"  ;
		std::cout.flush() ;

		if(temp_np <= _max_np)
		{
			++temp_np ;
		}
		if(temp_nq <= _max_nq)
		{
			++temp_nq ;
		}
	}

	engClose(ep); 
	return false ;
}

void rational_fitter_matlab::set_parameters(const arguments& args)
{
	_max_np = args.get_float("np", 10) ;
	_max_nq = args.get_float("nq", 10) ;
	_min_np = args.get_float("min-np", _max_np) ;
	_min_nq = args.get_float("min-nq", _max_nq) ;	
}
		
bool rational_fitter_matlab::fit_data(const rational_data* d, int np, int nq, rational_function* r) 
{

	// Multidimensional coefficients
	std::vector<double> Pn ; Pn.reserve(d->dimY()*np) ;
	std::vector<double> Qn ; Qn.reserve(d->dimY()*nq) ;

	for(int j=0; j<d->dimY(); ++j)
	{
		if(!fit_data(d, np, nq, j, r))
			return false ;

		for(int i=0; i<np; ++i) { Pn.push_back(r->getP(i)) ; }
		for(int i=0; i<nq; ++i) { Qn.push_back(r->getQ(i)) ; }
	}

	r->update(Pn, Qn) ;
	return true ;
}

// dat is the data object, it contains all the points to fit
// np and nq are the degree of the RP to fit to the data
// y is the dimension to fit on the y-data (e.g. R, G or B for RGB signals)
// the function return a ration BRDF function and a boolean
bool rational_fitter_matlab::fit_data(const rational_data* d, int np, int nq, int ny, rational_function* r) 
{
	// Size of the problem
	int N = np+nq ;
	int M = d->size() ;
	
	// Get the maximum value in data to scale the input parameter space
	// so that it reduces the values of the polynomial
	vec dmax = d->max() ;
	
	// Matrices of the problem
	Eigen::MatrixXd G (N, N) ;
	Eigen::VectorXd g (N) ;
	Eigen::MatrixXd CI(2*M, N) ;
	Eigen::VectorXd ci(2*M) ;

	// Select the size of the result vector to
	// be equal to the dimension of p + q
	for(int j=0; j<N; ++j)
	{
		for(int i=0; i<N; ++i)
		{
			if(i == j)
				G(i, j) = 1.0 ; 
			else
				G(i,j) = 0.0;
		}

		g(j) = 0.0 ;
	}
	
	// Each constraint (fitting interval or point
	// add another dimension to the constraint
	// matrix
	for(int i=0; i<M; ++i)	
	{		
		// Norm of the row vector
		double a0_norm = 0.0 ;
		double a1_norm = 0.0 ;

		vec xi = d->get(i) ;
		for(int k=0; k<d->dimX(); ++k)
		{
			xi[k] /= dmax[k] ;
		}

		// A row of the constraint matrix has this 
		// form: [p_{0}(x_i), .., p_{np}(x_i), -f(x_i) q_{0}(x_i), .., -f(x_i) q_{nq}(x_i)]
		// For the lower constraint and negated for 
		// the upper constraint
		for(int j=0; j<N; ++j)
		{
			// Filling the p part
			if(j<np)
			{
				const double pi = r->p(xi, j) ;
				a0_norm += pi*pi ;
				a1_norm += pi*pi ;
				// Updating Eigen matrix
				CI(2*i+0, j) =  pi ;
				CI(2*i+1, j) = -pi ;
			}
			// Filling the q part
			else
			{
				vec yl, yu ; 
				d->get(i, yl, yu) ;

				const double qi = r->q(xi, j-np) ;
				a0_norm += qi*qi * (yu[ny]*yu[ny]) ;
				a1_norm += qi*qi * (yl[ny]*yl[ny]) ;
				
				// Updating Eigen matrix
				CI(2*i+0, j) = -yu[ny] * qi ;
				CI(2*i+1, j) =  yl[ny] * qi ;
			}
		}
	
		// Set the c vector, will later be updated using the
		// delta parameter.
		ci(2*i+0) = -sqrt(a0_norm) ;
		ci(2*i+1) = -sqrt(a1_norm) ;
	}
	
	// Update the ci column with the delta parameter
	// (See Celis et al. 2007 p.12)
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(CI);
	const double sigma_m = svd.singularValues()(std::min(2*M, N)-1) ;
	const double sigma_M = svd.singularValues()(0) ;
	double delta = sigma_m / sigma_M ;
	
	if(std::isnan(delta) || (std::abs(delta) == std::numeric_limits<double>::infinity()))
	{
#ifdef DEBUG
		std::cerr << "<<ERROR>> delta factor is NaN of Inf" << std::endl ;
#endif
		return false ;
	}
	else if(delta == 0.0)
	{
		delta = 1.0 ;
	}

#ifdef DEBUG
	std::cout << "<<DEBUG>> delta factor: " << sigma_m << " / " << sigma_M << " = " << delta << std::endl ;
#endif
	for(int i=0; i<2*M; ++i)	
	{		
		ci(i) = ci(i) * delta ; 
	}
#ifdef DEBUG
	std::cout << "CI = " << CI << std::endl << std::endl ;
#endif

	// Create the MATLAB defintion of objects
	// MATLAB defines a quad prog as
	//   1/2 x' H x + f' x with A x <= b 
	//
	mxArray *H, *f, *A, *b, *x, *flag;
	H = mxCreateDoubleMatrix(  N, N, mxREAL);
	f = mxCreateDoubleMatrix(  N, 1, mxREAL);
	A = mxCreateDoubleMatrix(2*M, N, mxREAL);
	b = mxCreateDoubleMatrix(2*M, 1, mxREAL);

	memcpy((void *)mxGetPr(H), (void *) G.data(),  N*N*sizeof(double));
	memcpy((void *)mxGetPr(f), (void *) g.data(),  N*sizeof(double));
	memcpy((void *)mxGetPr(A), (void *)CI.data(), (2*M)*N*sizeof(double));
	memcpy((void *)mxGetPr(b), (void *)ci.data(), (2*M)*sizeof(double));

	engPutVariable(ep, "H", H);
	engPutVariable(ep, "f", f);
	engPutVariable(ep, "A", A);
	engPutVariable(ep, "b", b);
	
	char* output = new char[BUFFER_SIZE+1];
	output[BUFFER_SIZE] = '\0';
	engOutputBuffer(ep, output, BUFFER_SIZE) ;
#ifdef DEBUG
	engEvalString(ep, "display(H)");
	std::cout << output << std::endl ;
	engEvalString(ep, "display(f)");
	std::cout << output << std::endl ;
	engEvalString(ep, "display(A)");
	std::cout << output << std::endl ;
	engEvalString(ep, "display(b)");
	std::cout << output << std::endl ;
#endif
	engEvalString(ep, "[x, fval, flag] = quadprog(H,f,A,b);");
#ifdef DEBUG
	std::cout << output << std::endl ;
	engEvalString(ep, "display(x)");
	std::cout << output << std::endl ;
	engEvalString(ep, "display(flag)");
#endif
	mxDestroyArray(H);
	mxDestroyArray(f);
	mxDestroyArray(A);
	mxDestroyArray(b);

	x    = engGetVariable(ep, "x") ;
	flag = engGetVariable(ep, "flag") ;
	if(x != NULL)
	{
		if(flag != NULL)
		{
			if(mxGetScalar(flag) != 1)
			{
				mxDestroyArray(x);
				mxDestroyArray(flag);
			
#ifdef DEBUG
				std::cerr << "<<ERROR>> flag is not equal to 1" << std::endl ;
#endif
				return false ;
			}

			double  total = 0.0;
			double* val = (double*)mxGetData(x) ;
			std::vector<double> a, b;
			for(int i=0; i<N; ++i)
			{
				total += val[i]*val[i] ;
				if(i < np)
				{
					a.push_back(val[i] / r->p(dmax, i)) ;
				}
				else
				{
					b.push_back(val[i] / r->q(dmax, i-np)) ;
				}
			}
			r->update(a, b) ;

			mxDestroyArray(x);
			mxDestroyArray(flag);
			return total > 0.0 ;
		}
		else
		{
#ifdef DEBUG
			std::cerr << "<<ERROR>> unable to gather result flag" << std::endl ;
#endif
			return false ;
		}
	}
	else 
	{
#ifdef DEBUG
		std::cerr << "<<ERROR>> unable to gather result x" << std::endl ;
#endif
		return false ;
	}
}

Q_EXPORT_PLUGIN2(rational_fitter_matlab, rational_fitter_matlab)
