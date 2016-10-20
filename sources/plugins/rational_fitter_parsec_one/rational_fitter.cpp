#include <Eigen/SVD>
#include <QuadProg++.hh>
#include <core/plugins_manager.h>

#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <cmath>

#include <dplasma.h>
#include <data_dist/matrix/two_dim_rectangle_cyclic.h>
#include <gesvdm.h>

#include "rational_fitter.h"

#ifdef WIN32
#define isnan(X) ((X != X))
#endif

#include <core/common.h>

ALTA_DLL_EXPORT fitter* provide_fitter()
{
    return new rational_fitter_parsec_multi();
}

rational_fitter_parsec_multi::rational_fitter_parsec_multi() : _boundary(1.0), _args(NULL)
{
}

rational_fitter_parsec_multi::~rational_fitter_parsec_multi()
{
    dague_fini(&_dague);
}

void rational_fitter_parsec_multi::set_parameters(const arguments& args)
{
    _max_np = args.get_float("np", 10);
    _max_nq = args.get_float("nq", _max_np);
    _min_np = args.get_float("min-np", _max_np);
    _min_nq = args.get_float("min-nq", _min_np);
 
    _max_np = std::max<int>(_max_np, _min_np);
    _max_nq = std::max<int>(_max_nq, _min_nq);
    
    _boundary = args.get_float("boundary-constraint", 1.0f);
    _nbcores = args.get_int( "nbcores", 2 );
    _args = &args;

    {
	int argc = 1;
	char **argv = (char**)malloc((argc+1)*sizeof(char*));

	argv[0] = strdup( "./manao" );
	argv[1] = NULL;

	_dague = dague_init(_nbcores, &argc, &argv);

	free(argv[0]);
	free(argv);
    }
}


bool rational_fitter_parsec_multi::fit_data(const ptr<data>& dat, ptr<function>& fit, const arguments &args)
{
    std::cerr << "Entering first level of fit_data" << std::endl;

    ptr<rational_function> r = dynamic_pointer_cast<rational_function>(fit);
    const ptr<vertical_segment>& d = dynamic_pointer_cast<vertical_segment>(dat);
    if(!r || !d
       || d->confidence_interval_kind() != vertical_segment::ASYMMETRICAL_CONFIDENCE_INTERVAL)
	{
	    std::cerr << "<<ERROR>> not passing the correct class to the fitter" << std::endl;
	    return false;
	}

    // I need to set the dimension of the resulting function to be equal
    // to the dimension of my fitting problem

    // This doesn't work for dimension greater than 1 for now.
    //assert( d->dimY() == 1 );

    r->setDimX(d->dimX());
    r->setDimY(d->dimY());
    r->setMin(d->min());
    r->setMax(d->max());

    std::cout << "<<INFO>> np in  [" << _min_np << ", " << _max_np
	      << "] & nq in [" << _min_nq << ", " << _max_nq << "]" << std::endl;

    const int step = args.get_int("np-step", 1);
    int i;

    for(i=std::max<int>(2,_min_np); i<=_max_np; i+=step )
	{
	    // double min_delta  = std::numeric_limits<double>::max();
	    // double min_l2_dist= std::numeric_limits<double>::max();
	    // double mean_delta = 0.0;
	    // int nb_sol_found  = 0;
	    // int nb_sol_tested = 0;
	    timer time;
	    int np;

	    std::cout << "<<INFO>> fit using np+nq = " << i << std::endl ;
	    std::cout.flush() ;

	    time.start();

	    // i=np+nq => i-1 independant pb
	    std::cout << r.getcounter() << std::endl;
	    if(fit_data(d.get(), i-1, r.get(), np))
		{
		    time.stop();

		    std::cout << r.getcounter() << std::endl;
		    std::cout << "<<INFO>> got a fit using np = " << np << " & nq =  " << i-np << "      " << std::endl;
		    std::cout << "<<INFO>> " << *(r.get()) << std::endl;
		    std::cout << "<<INFO>> it took " << time << std::endl;

		    return true;
		}
	}

    std::cerr << "Exiting first level of fit_data" << std::endl;
    return false;
}

bool rational_fitter_parsec_multi::fit_data(vertical_segment *d, int N, rational_function *r, int &np)
{
    std::cerr << "Entering Second level of fit_data for (np+nq)=" << N+1 << std::endl;

    // Size of the problem
    const int M = d->size();
    const int twoM = 2 * M;
    bool rc = false;

    /*
     * Create tree for QR reduction
     */
    {
	gesvdm2_args_t args;
	tiled_matrix_desc_t *panel;
	dplasma_qrtree_t qrtree;
	two_dim_block_cyclic_t *matrixCI = (two_dim_block_cyclic_t*)malloc(sizeof(two_dim_block_cyclic_t));
	two_dim_block_cyclic_t *vectorCI = (two_dim_block_cyclic_t*)malloc(sizeof(two_dim_block_cyclic_t));
	int nb = N + 1;
	int mb = nb; //twoM; //nb;
	//int mb = std::max<int>(nb, 400);

	//	assert( (mb == nb) && (mb%2 == 0) );
	if (mb%2 == 1)
	    mb++;

	memset(matrixCI, 0, sizeof(two_dim_block_cyclic_t));
	memset(vectorCI, 0, sizeof(two_dim_block_cyclic_t));

	two_dim_block_cyclic_init( matrixCI, matrix_RealDouble, matrix_Lapack,
				   1, 0, mb, nb, twoM, nb * N, 0, 0,
				   twoM, nb * N, 1, 1, 1);

	/* Overwrite llm top make sure we use the correct one and not a multiple of N */
	matrixCI->super.llm = twoM;
	//matrixCI->mat = (double*)malloc( nb*N * twoM * sizeof(double));
	matrixCI->mat = NULL;

	panel = tiled_matrix_submatrix( (tiled_matrix_desc_t *)matrixCI,
					0, 0, twoM, (twoM > nb ? nb : twoM) );

	std::cerr << "<<DEBUG>> M: " << M << " N: "<< N << std::endl;

	dplasma_hqr_init( &qrtree, PlasmaNoTrans, panel,
			  DPLASMA_FLAT_TREE, DPLASMA_FIBONACCI_TREE,
			  4, 1, 0, 0 );

	r->setSize( nb, nb );
	for(int j=0; j<d->dimY(); j++)
	    r->get(j)->resize( nb, nb );

	args.qrtree = &qrtree;
	args.ib = ((32 < nb) ? 32 : nb);
	args.M = twoM;
	args.n = N;
	args.ndim = d->dimY();
	args.A  = (tiled_matrix_desc_t *)matrixCI;
	args.dataptr = d;
	args.rfpm    = this;
	args.rfptr   = r;

	args.fillp = &fill_p;
	args.fillq = &fill_q;
	args.solve_init     = &solve_init;
	args.solve          = &solve_wrapper;
	args.solve_finalize = &solve_finalize;

	args.pbs   = (subproblem_t *)calloc(N, sizeof(subproblem_t));

	std::cerr << "<<DEBUG>> Start gesvd" << std::endl;
	dplasma_dgesvdm2( _dague, &args );

	std::cerr << "<<DEBUG>> End of gesvd" << std::endl;
	free(panel);
	dplasma_hqr_finalize( &qrtree );

	dague_ddesc_destroy((dague_ddesc_t*)matrixCI);
	dague_ddesc_destroy((dague_ddesc_t*)vectorCI);
	std::cerr << "<<DEBUG>> Data freed" << std::endl;

	free(matrixCI);
	free(vectorCI);

	{
	    subproblem_t *spb = args.pbs;
	    int nb_sol_found  = 0;
	    int min_sol = -1;
	    double min_delta   = std::numeric_limits<double>::max();
	    double min_l2_dist = std::numeric_limits<double>::max();
	    double mean_delta = 0.0;

	    for (int i=0; i<N; i++, spb++) {
		
		if ( spb->isfitted == 1 ) {
		    rc = true;
		    nb_sol_found ++;
		    mean_delta += spb->delta[0];
		    
		    std::cout << "<<INFO>> found a solution with np=" << (i+1)
			      << ", nq= " << (N-i)
			      << ", delta= " << spb->delta[0]
			      << ", l2_dist= " << spb->l2_dist
			      << std::endl;
		    
		    if(spb->delta[0] < min_delta)
			{
			    min_delta   = spb->delta[0] ;
			    min_l2_dist = spb->l2_dist ;
			    min_sol = i;
			}
		}
	    }
	    std::cout << std::endl;
	
	    if (rc) {
		np = min_sol+1;
		r->setSize( min_sol+1, N-min_sol );
		r->update( (rational_function*)(args.pbs[min_sol].rfptr) );
		
		std::cout << "<<INFO>> mean delta = " << mean_delta/nb_sol_found << std::endl;
		std::cout << "<<INFO>>  min delta = " << min_delta << std::endl;
		std::cout << "<<INFO>>  min l2 dist = " << min_l2_dist << std::endl;
	    }
	}
    }

    std::cerr << "Exiting Second level of fit_data for (np+nq)=" << N+1 << std::endl;
    return rc;
}

/**
 * @param args[in] 
 * @param ny[in] 
 *         Specifies the output direction
 *
 * @param i0[in] 
 *         Specifies the firs tow to generate
 *
 * @param M[in] 
 *         Specifies the number of row to generate in the matrix P.
 *
 * @param P[out]
 *         Array of dimension (args->n+1) * M
 *         Matrix to store the polynomial evaluation of the coefficients
 *
 * @param ldp[in] 
 *         Leading dimension of the matrix P. ldp == args-n+1.
 *
 */
void rational_fitter_parsec_multi::fill_p(const gesvdm2_args_t *args, int ny, int i0, int M,
					  double *P, int ldp)
{
    const vertical_segment *d    = (vertical_segment *)(args->dataptr);
    rational_function_1d   *rf1d = ((rational_function *)(args->rfptr))->get(ny);
    double *lP = P;
    vec xi;
    const int N = args->n + 1;
    int i, j;

    /* Check row major format with exact number of columns */
    assert(ldp == N);

    for( i=0; i<M; i++, i0++) {
	xi = d->get(i0);

	// A row of the P matrix has this
	// form: [p_{0}(x_i), .., p_{np}(x_i)]
	// The lower and upper constraint are later computed
	for(j=0; j<N; j++, lP++) {
	    *lP = rf1d->p(xi, j);
	}
    }
}

void rational_fitter_parsec_multi::fill_q(const gesvdm2_args_t *args, int ny, int i0, int M,
					  double *Q, int ldq)
{
    const vertical_segment *d    = (vertical_segment *)(args->dataptr);
    rational_function_1d   *rf1d = ((rational_function *)(args->rfptr))->get(ny);
    double *lQ = Q;
    vec xi, yl, yu;
    const int N = args->n + 1;
    int i, j;

    /* Check row major format with exact number of columns */
    assert( ldq == (N + 2) );

    for( i=0; i<M; i++, i0++) {
	d->get(i0, xi, yl, yu);

	*lQ = yu[ny]; lQ++;
	*lQ = yl[ny]; lQ++;

	// A row of the Q matrix has this
	// form: [q_{0}(x_i), .., q_{nq}(x_i)]
	// The lower and upper constraint are later computed
	for(j=0; j<N; j++, lQ++) {
	    *lQ = rf1d->q(xi, j);
	}
    }
}

int rational_fitter_parsec_multi::test_all_constraint( const vertical_segment     *data,
						       const rational_function_1d *r,
						       int ny )
{
    int n;
    for(n=0; n<data->size(); ++n)
    {
	vec x, yl, yu;
	data->get(n, x, yl, yu);

	vec y = r->value(x);
	if( (y[0] < yl[ny]) ||
	    (y[0] > yu[ny]) )
	{
	    return 0;
	}
    }
    return 1;
}

/**
 * Take as input the full constraint matrix CI for a dimension ny, and
 * the associated ci vector with the 2-norm of each row scaled by delta.
 * Apply the quadratic program and test against all the constraints.
 * Return true if all constraints match, false otherwise.
 */
void rational_fitter_parsec_multi::solve_init( const gesvdm2_args_t *args, subproblem_t *spb )
{
    rational_fitter_parsec_multi *rfpm = (rational_fitter_parsec_multi*)(args->rfpm) ;
    rational_function            *rf   = (rational_function*)(args->rfptr);
    rational_function            *rf2  = NULL;
    int np = spb->np;
    int nq = spb->nq;

    rf2 = dynamic_cast<rational_function*>(plugins_manager::get_function(*(rfpm->_args)));

    std::cout << rf2 << ": Start initPbRf" << std::endl;
    if(!rf)
	{
	    std::cerr << "<<ERROR>> unable to obtain a rational function from the plugins manager" << std::endl;
	    throw;
	}
    rf2->setParametrization(rf->input_parametrization());
    rf2->setParametrization(rf->output_parametrization());
    rf2->setDimX(rf->dimX()) ;
    rf2->setDimY(rf->dimY()) ;
    rf2->setMin(rf->min()) ;
    rf2->setMax(rf->max()) ;
	
    // Set the rational function size
    rf2->setSize(np, nq);

    spb->rfptr = rf2;
}


/**
 * Take as input the full constraint matrix CI for a dimension ny, and
 * the associated ci vector with the 2-norm of each row scaled by delta.
 * Apply the quadratic program and test against all the constraints.
 * Return true if all constraints match, false otherwise.
 */
int rational_fitter_parsec_multi::solve_wrapper( const gesvdm2_args_t *args, subproblem_t *pb, int M, int ny,
						 double *CIptr, double *ciptr )
{
    const vertical_segment *d    = (const vertical_segment*)(args->dataptr);
    rational_function      *rf   = (rational_function*)(pb->rfptr);
    rational_function_1d   *rf1d = rf->get(ny);
    const int np = pb->np;
    const int nq = pb->nq;
    const int N  = np + nq;

    // Compute the solution
    Eigen::MatrixXd CE(N, 0);
    Eigen::VectorXd ce(0);
    Eigen::MatrixXd G (N, N); G.setIdentity();
    Eigen::VectorXd g (N); g.setZero();
    Eigen::VectorXd x (0.0, N);

    Eigen::MatrixXd CI = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>::Map(CIptr, N, M);
    Eigen::Map<Eigen::VectorXd> ci(ciptr, M);

    // Select the size of the result vector to
    // be equal to the dimension of p + q
    double cost = QuadProgPP::solve_quadprog(G, g, CE, ce, CI, ci, x);
    bool solves_qp = !(cost == std::numeric_limits<double>::infinity());

    if(solves_qp) {
	std::cout << "<<INFO>> got solution for pb with np=" << pb->np << ", nq=" << pb->nq << std::endl;

	// Recopy the vector d
	vec p(np), q(nq);
	double norm = 0.0;

	for(int i=0; (i<N) & solves_qp; ++i) {
	    const double v = x[i];

	    solves_qp = solves_qp && !isnan(v)
		&& (v != std::numeric_limits<double>::infinity());

	    norm += v*v;
	    if(i < np) {
		p[i] = v;
	    }
	    else {
		q[i - np] = v;
	    }
	}

	if (solves_qp) {
	    std::cout << "<<INFO>> got solution to second step for pb with np=" << pb->np << ", nq=" << pb->nq << std::endl;

	    // Rq: doesn't need protection in // since it should be working on independant vectors
	    rf1d->update(p, q);
	    solves_qp = (test_all_constraint( d, rf1d, ny ) == 1);
	}
    }
    else {
	std::cerr << "<<DEBUG>> Didn't get solution to the pb with np=" << pb->np << ", nq=" << pb->nq << std::endl;
    }

    return solves_qp;
}

int rational_fitter_parsec_multi::solve_finalize( const gesvdm2_args_t *args, subproblem_t *pb )
{
    const vertical_segment *d  = (const vertical_segment*)(args->dataptr);
    rational_function      *rf = (rational_function*)(pb->rfptr);

    static ptr<data> dat((data *)d);
    
    std::cout << "d=" << d << std::endl;
    std::cout << "<<INFO>> got solution to final step for pb with np=" << pb->np << ", nq=" << pb->nq << std::endl;
    pb->linf_dist = rf->Linf_distance(dat);
    pb->l2_dist   = rf->L2_distance(dat);
    std::cout << "<<INFO>> got solution with np=" << pb->np << ", nq=" << pb->nq
	      << " Linf: " << pb->linf_dist << " L2: " << pb->l2_dist << std::endl;
    std::cout << "d=" << d << std::endl;

    return pb->isfitted;
}
