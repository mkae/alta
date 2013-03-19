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

bool rational_fitter_dca::fit_data(const data* dat, function* fit)
{
	rational_function* r = dynamic_cast<rational_function*>(fit) ;
	const vertical_segment* d = dynamic_cast<const vertical_segment*>(dat) ;
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
	r->setMin(d->min()) ;
	r->setMax(d->max()) ;

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

void rational_fitter_dca::set_parameters(const arguments& args)
{
	_max_np = args.get_float("np", 10) ;
	_max_nq = args.get_float("nq", 10) ;
	_min_np = args.get_float("min-np", _max_np) ;
	_min_nq = args.get_float("min-nq", _max_nq) ;	
}
		
bool rational_fitter_dca::fit_data(const data* d, int np, int nq, rational_function* r)
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


// Bootstrap the DCA algorithm with the Papamarkos fitting
// algorithm [Papamarkos 1988]
// \todo Finish the Papamarkos implementation
void bootstrap(const data* d, int np, int nq, rational_function* fit, double& delta)
{
}

// dat is the data object, it contains all the points to fit
// np and nq are the degree of the RP to fit to the data
// y is the dimension to fit on the y-data (e.g. R, G or B for RGB signals)
// the function return a ration BRDF function and a boolean
bool rational_fitter_dca::fit_data(const data* d, int np, int nq, int ny, rational_function* r)
{
	// Size of the problem
	int N  = np+nq+1 ;
	int M  = d->size() ;
	int nY = d->dimY();	

    // Bootstrap the delta and rational function using the Papamarkos
    // algorithm.
    double delta = 0.0;
    bootstrap(d, np, nq, r, delta);

	// Create the MATLAB defintion of objects
	// MATLAB defines a linear prog as
	//   min f' x with A x <= b 
	//
	mxArray *f, *A, *b, *x, *flag;
	f = mxCreateDoubleMatrix(N*nY,    1, mxREAL);
	A = mxCreateDoubleMatrix( 2*M, nY*N, mxREAL);
	b = mxCreateDoubleMatrix( 2*M,    1, mxREAL);

	engPutVariable(ep, "f", f);
	engPutVariable(ep, "A", A);
	engPutVariable(ep, "b", b);

	// Matrices of the problem in Eigen format
	Eigen::VectorXd g (nY*N) ;
	Eigen::MatrixXd CI(2*M, nY*N) ;
	Eigen::VectorXd ci(2*M) ;

    double delta_k = delta;

	// Loop until you get a converge solution \delta > \delta_k
	// \todo add the correct looping condition
	while(true)
	{
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
			// Norm of the row vector
			double a0_norm = 0.0 ;
			double a1_norm = 0.0 ;

			vec xi = d->get(i) ;

			// For each input data i \in M, the following constraints have to
			// be fulfilled:
			//   [ f_i + \delta_k] * q_i - p_i + qk_i \delta >= 0
			//   [-f_i + \delta_k] * q_i + p_i + qk_i \delta >= 0
			//                                            pi >= 0
			//                                            qi >= 0
			//
			// The constraints matrix has the following form
			// [-p_j(x_i) ..-p_j(x_i), [ f_i+\delta_k]*q_j(x_i) ..[ fi+\delta_k]*q_j(x_i), qk(x_i)]
			// [ p_j(x_i) .. p_j(x_i), [-f_i+\delta_k]*q_j(x_i) ..[-fi+\delta_k]*q_j(x_i), qk(x_i)]
			//
			for(int j=0; j<N-1; ++j)
			{
				// Filling the p part
				if(j<np)
				{
					const double pi = r->p(xi, j) ;

					// Updating Eigen matrix
					for(int y=0; y<nY; ++y)
					{
						CI(2*(nY*i + y)+0, j) = -pi ;
						CI(2*(nY*i + y)+1, j) =  pi ;
					}
				}
				// Filling the q part
				else if(j<np+nq)
				{
                    vec value = d->get(i) ;

					const double qi = r->q(xi, j-np) ;

					// Updating Eigen matrix
					for(int y=0; y<nY; ++y)
					{
                        CI(2*(nY*i + y)+0, j) = (delta_k+value[y]) * qi ;
                        CI(2*(nY*i + y)+1, j) = (delta_k-value[y]) * qi ;
					}
				}
				else
				{
					// Last column of the constraint matrix
                    vec qk = r->q(xi) ;
					for(int y=0; y<nY; ++y)
					{
						CI(2*(nY*i + y)+0, j) = qk[y] ;
						CI(2*(nY*i + y)+1, j) = qk[y] ;
					}
				}
			}


			// Set the c vector
			for(int y=0; y<nY; ++y)
			{
				ci(2*(nY*i+y)+0) = 0.0 ;
				ci(2*(nY*i+y)+1) = 0.0 ;
			}
		}


		// Copy the data to matlab and execute the linear program
		//
		memcpy((void *)mxGetPr(f), (void *) g.data(),  N*sizeof(double));
		memcpy((void *)mxGetPr(A), (void *)CI.data(), (2*M)*N*sizeof(double));
		memcpy((void *)mxGetPr(b), (void *)ci.data(), (2*M)*sizeof(double));

		char* output = new char[BUFFER_SIZE+1];
		output[BUFFER_SIZE] = '\0';
		engOutputBuffer(ep, output, BUFFER_SIZE) ;
#ifdef DEBUG
		engEvalString(ep, "display(f)");
		std::cout << output << std::endl ;
		engEvalString(ep, "display(A)");
		std::cout << output << std::endl ;
		engEvalString(ep, "display(b)");
		std::cout << output << std::endl ;
#endif

		engEvalString(ep, "[x, fval, flag] = linprog(f,A,b);");
#ifdef DEBUG
		std::cout << output << std::endl ;
#endif

		x    = engGetVariable(ep, "x") ;
		flag = engGetVariable(ep, "flag") ;

		// \todo remove when correct looping condition ready
		break;
	}

	mxDestroyArray(f);
	mxDestroyArray(A);
	mxDestroyArray(b);
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
					a.push_back(val[i]) ;
				}
				else
				{
					b.push_back(val[i]) ;
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

Q_EXPORT_PLUGIN2(rational_fitter_dca, rational_fitter_dca)