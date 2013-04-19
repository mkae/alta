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

fitter* provide_fitter()
{
    return new rational_fitter_leastsquare();
}

data* rational_fitter_leastsquare::provide_data() const
{
	return new vertical_segment() ;
}

function* rational_fitter_leastsquare::provide_function() const 
{
	return new rational_function() ;
}

rational_fitter_leastsquare::rational_fitter_leastsquare() : QObject()
{
}
rational_fitter_leastsquare::~rational_fitter_leastsquare() 
{
}

bool rational_fitter_leastsquare::fit_data(const data* dat, function* fit, const arguments &args)
{
	rational_function* r = dynamic_cast<rational_function*>(fit) ;
	const vertical_segment* d = dynamic_cast<const vertical_segment*>(dat) ;
	if(r == NULL || d == NULL)
	{
		std::cerr << "<<ERROR>> not passing the correct class to the fitter" << std::endl ;
		return false ;
	}

	// I need to set the dimension of the resulting function to be equal
	// to the dimension of my fitting problem
	r->setDimX(d->dimX()) ;
	r->setDimY(d->dimY()) ;
	r->setMin(d->min()) ;
	r->setMax(d->max()) ;

	std::cout << "<<INFO>> np =" << _np << "& nq =" << _nq  << std::endl ;


	QTime time ;
	time.start() ;


	if(fit_data(d, _np, _nq, r))
	{
		int msec = time.elapsed() ;
		int sec  = (msec / 1000) % 60 ;
		int min  = (msec / 60000) % 60 ;
		int hour = (msec / 3600000) ;
		std::cout << "<<INFO>> got a fit" << std::endl ;
		std::cout << "<<INFO>> it took " << hour << "h " << min << "m " << sec << "s" << std::endl ;

		return true ;
	}

	std::cout << "<<INFO>> fit failed\r"  ;
	std::cout.flush() ;

	return false ;
}

void rational_fitter_leastsquare::set_parameters(const arguments& args)
{
	_np = args.get_float("np", 10) ;
	_nq = args.get_float("nq", 10) ;
  _max_iter = args.get_float("max-iter", 1) ;
}
		
bool rational_fitter_leastsquare::fit_data(const vertical_segment* d, int np, int nq, rational_function* r) 
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
bool rational_fitter_leastsquare::fit_data(const vertical_segment* d, int np, int nq, int ny, rational_function* r) 
{
  using namespace Eigen;
  
  double scale = 1;
  
  MatrixXd D(d->size(), np+nq);
  VectorXd Y(d->size());
  for(int i=0; i<d->size(); ++i) 
  {
    vec yl, yu ;
    d->get(i, yl, yu);
    const double y  = 0.5*(yl[ny] + yu[ny]) ;
    Y[i] = y;
    vec x = d->get(i);
    VectorXd::Map(&x[0], 1) /= scale;
    // A row of the constraint matrix has this 
    // form: [p_{0}(x_i), .., p_{np}(x_i), q_{0}(x_i), .., q_{nq}(x_i)]
    for(int j=0; j<np+nq; ++j)
    {
      if(j<np)
      {
        const double pi = r->p(x, j) ;
        D(i, j) =  pi ;
      }
      else
      {
        const double qi = r->q(x, j-np) ;
        D(i,j) = qi ;
      }
    }
  }

  VectorXd pq(np+nq);
  
  // initialize with a constant numerator:
  pq.head(np).setZero();
  pq(0) = 1;
  
  // Alternate between fixed numerator and fixed denominator
  // By default we iterate only once (usually successive iteration introduced more errors)
  for(int iter=0; iter<_max_iter; ++iter)
  {
    
    // Step 1 fit 1/f_i using 1/q only
    {
      // Theoretically best weighting scheme to approcimate the true LS problem, which is: sum_i (1/q(x_i) - f_i)^2
      MatrixXd A = Y.cwiseAbs2().asDiagonal() * D.rightCols(nq);
      VectorXd b = Y.asDiagonal() * (D.leftCols(np) * pq.head(np));
      VectorXd x = A.colPivHouseholderQr().solve(b);
      
      // Sometimes, this scheme works better:
//       MatrixXd A = Y.asDiagonal() * D.rightCols(nq);
//       VectorXd b = D.leftCols(np) * pq.head(np);
//       VectorXd x = A.colPivHouseholderQr().solve(b);
      
      pq.tail(nq) = x;
      
      // Statistics on the problem:
//       VectorXd sv = A.jacobiSvd().singularValues();
//       std::cout << "Singular values: " << sv.transpose() << std::endl;
//       std::cout << "Rcond: " << sv(0)/sv(nq-1) << std::endl;
//       std::cout << "<<INFO>> LS error " << (A*x-b).norm()/b.norm() << std::endl;
    }
    
    VectorXd res = ((D.leftCols(np) * pq.head(np)).array() / (D.rightCols(nq) * pq.tail(nq)).array() - Y.array());
/*
  std::cout << "<<INFO>> Real LS norm (1): "
		 << res.norm() / Y.norm()
		 << " ; inf " << res.lpNorm<Infinity>() / Y.lpNorm<Infinity>()
		 << " ; L1  " << res.lpNorm<1>() / Y.lpNorm<1>() << std::endl;
*/
    // Step 2 fit f_i using p/q with q fix
    {
      // This scheme gives more weights to small values: (not good)
  //     VectorXd V = D.rightCols(nq) * pq.tail(nq);
  //     MatrixXd A = Y.asDiagonal().inverse() * (V.asDiagonal().inverse() * D.leftCols(np));
  //     VectorXd b = VectorXd::Ones(Y.size());
  //     VectorXd x = A.colPivHouseholderQr().solve(b);
      
      VectorXd Q = D.rightCols(nq) * pq.tail(nq);
      MatrixXd A = Q.asDiagonal().inverse() * D.leftCols(np);
      VectorXd x = A.colPivHouseholderQr().solve(Y);
      
      pq.head(np) = x;

      // Statistics on the problem:
//       VectorXd sv = A.jacobiSvd().singularValues();
//       std::cout << "Singular values: " << sv.transpose() << std::endl;
//       std::cout << "Rcond: " << sv(0)/sv(np-1) << std::endl;
//       std::cout << "<<INFO>> LS error " << (A*x-Y).norm()/Y.norm() << std::endl;
    }
  } // iterations
  
  std::vector<double> p(np), q(nq);
  Eigen::VectorXd::Map(&p[0], np) = pq.head(np);
  Eigen::VectorXd::Map(&q[0], nq) = pq.tail(nq);
  
  double s = 1;
  for(int i=0; i<np; ++i)
  {
    p[i] /= s;
    s *= scale;
  }
  s = 1;
  for(int i=0; i<nq; ++i)
  {
    q[i] /= s;
    s *= scale;
  }
  
  // Evaluate true LS error: sum_i (p(x_i)/q(x_i) - f_i)^2
  VectorXd res = ((D.leftCols(np) * pq.head(np)).array() / (D.rightCols(nq) * pq.tail(nq)).array() - Y.array());
  std::cout << "<<INFO>> L_2 "
    << res.norm() / Y.norm()
    << " ; L_inf " << res.lpNorm<Infinity>() / Y.lpNorm<Infinity>()
    << " ; L_1  " << res.lpNorm<1>() / Y.lpNorm<1>() << std::endl;

  r->update(p, q) ;
  std::cout << "<<INFO>> got solution " << *r << std::endl ;
  return true;

}

Q_EXPORT_PLUGIN2(rational_fitter_leastsquare, rational_fitter_leastsquare)
