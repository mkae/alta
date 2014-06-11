#include <Eigen/SVD>
#include <Array.hh>
#include <QuadProg++.hh>

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
    return new rational_fitter_parsec_one();
}

rational_fitter_parsec_one::rational_fitter_parsec_one() : _boundary(1.0)
{
    //args.get_int( "nbcores" );

    nbcores = 2;
    int argc = 3;
    char **argv = (char**)malloc(argc*sizeof(char*));

    argv[0] = strdup( "./dague" );
    argv[1] = strdup( "--" );
    argv[2] = strdup( "--dague_dot" );

    dague = dague_init(nbcores, &argc, &argv);
    // {
    //     int p, nb_total_comp_threads = 0;
    //     for(p = 0; p < dague->nb_vp; p++) {
    //         nb_total_comp_threads += dague->virtual_processes[p]->nb_cores;
    //     }
    //     nbcores = nb_total_comp_threads;
    // }
}

rational_fitter_parsec_one::~rational_fitter_parsec_one()
{
    dague_fini(&dague);
}

bool rational_fitter_parsec_one::fit_data(const ptr<data>& dat, ptr<function>& fit, const arguments &args)
{
    ptr<rational_function> r = dynamic_pointer_cast<rational_function>(fit) ;
    const ptr<vertical_segment> d = dynamic_pointer_cast<vertical_segment>(dat) ;
    if(!r || !d)
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

    std::cout << "<<INFO>> np in  [" << _min_np << ", " << _max_np
              << "] & nq in [" << _min_nq << ", " << _max_nq << "]" << std::endl ;

    int temp_np = _min_np, temp_nq = _min_nq ;
    while(temp_np <= _max_np || temp_nq <= _max_nq)
        {
            timer time ;
            time.start() ;

            r->setSize(temp_np, temp_nq);
            if(fit_data(d, temp_np, temp_nq, r))
                {
                    time.stop() ;

                    std::cout << "<<INFO>> got a fit using np = " << temp_np << " & nq =  " << temp_nq << "      " << std::endl ;
                    std::cout << "<<INFO>> it took " << time << std::endl ;

                    return true ;
                }


            std::cout << "<<INFO>> fit using np = " << temp_np << " & nq =  " << temp_nq << " failed" << std::endl  ;
            std::cout.flush() ;

            if(temp_np == _max_np && temp_nq == _max_nq)
                {
                    return false;
                }

            if(temp_np < _max_np)
                {
                    ++temp_np ;
                }
            if(temp_nq < _max_nq)
                {
                    ++temp_nq ;
                }
        }
    return false ;
}

void rational_fitter_parsec_one::set_parameters(const arguments& args)
{
    _max_np = args.get_float("np", 10) ;
    _max_nq = args.get_float("nq", 10) ;
    _min_np = args.get_float("min-np", _max_np) ;
    _min_nq = args.get_float("min-nq", _max_nq) ;

    _max_np = std::max<int>(_max_np, _min_np);
    _max_nq = std::max<int>(_max_nq, _min_nq);

    _boundary = args.get_float("boundary-constraint", 1.0f);
}


bool rational_fitter_parsec_one::fit_data(const ptr<vertical_segment>& d, int np, int nq, const ptr<rational_function>& r)
{
    // For each output dimension (color channel for BRDFs) perform
    // a separate fit on the y-1D rational function.
    for(int j=0; j<d->dimY(); ++j)
        {
            rational_function_1d* rs = r->get(j);
            rs->resize(np, nq);

            if(!fit_data(d, np, nq, j, rs))
                {
                    return false ;
                }
        }

    return true ;
}

// dat is the data object, it contains all the points to fit
// np and nq are the degree of the RP to fit to the data
// y is the dimension to fit on the y-data (e.g. R, G or B for RGB signals)
// the function return a ration BRDF function and a boolean
bool rational_fitter_parsec_one::fit_data(const ptr<vertical_segment>& d, int np, int nq, int ny, rational_function_1d* r)
{
    // Size of the problem
    const int N = np+nq ;
    const int M = d->size() ;
    const int twoM = 2 * M;
    bool rc = false;

    // Matrices of the problem
    QuadProgPP::Matrix<double> CI(0.0, N, twoM) ;
    QuadProgPP::Vector<double> ci(0.0, twoM) ;

    /*
     * Create tree for QR reduction
     */
    {
        gesvdm_args_t args;
        tiled_matrix_desc_t *panel;
        dplasma_qrtree_t qrtree;
        two_dim_block_cyclic_t *matrixCI = (two_dim_block_cyclic_t*)malloc(sizeof(two_dim_block_cyclic_t));
        two_dim_block_cyclic_t *vectorCI = (two_dim_block_cyclic_t*)malloc(sizeof(two_dim_block_cyclic_t));

        memset(matrixCI, 0, sizeof(two_dim_block_cyclic_t));
        memset(vectorCI, 0, sizeof(two_dim_block_cyclic_t));

        two_dim_block_cyclic_init( matrixCI, matrix_RealDouble, matrix_Lapack,
                                   1, nbcores, 0, N, N, twoM, N, 0, 0,
                                   twoM, N, 1, 1, 1);
        two_dim_block_cyclic_init( vectorCI, matrix_RealDouble, matrix_Tile,
                                   1, nbcores, 0, N, 1, twoM, 1, 0, 0,
                                   twoM, 1, 1, 1, 1);
        /* Overwrite llm top make sure we use the correct one and not a multiple of N */
        matrixCI->super.llm = twoM;
        matrixCI->mat = CI[0];
        vectorCI->super.llm = twoM;
        vectorCI->mat = &(ci[0]);

        panel = tiled_matrix_submatrix( (tiled_matrix_desc_t *)matrixCI,
                                        0, 0, twoM, (twoM > N ? N : twoM) );

        std::cerr << "<<DEBUG>> M: " << M << " N: "<< N << std::endl;

        dplasma_hqr_init( &qrtree, PlasmaNoTrans, panel,
                          DPLASMA_FLAT_TREE, DPLASMA_FIBONACCI_TREE,
                          4, 1, 0, 0 );

        args.qrtree = &qrtree;
        args.ib = ((32 < N) ? 32 : N);
        args.M = M;
        args.np = np;
        args.nq = nq;
        args.ny = ny;
        args.A  = (tiled_matrix_desc_t *)matrixCI;
        args.ci = (tiled_matrix_desc_t *)vectorCI;
        args.data = d.get();
        args.rptr = r;
        args.fillp = &fillci_p;
        args.fillq = &fillci_q;
        args.solve = &solve_wrapper;
        args.rc = 0;

        std::cerr << "<<DEBUG>> Start gesvd" << std::endl;
        dplasma_dgesvd_multiple( dague, &args );
        rc = args.rc;

        std::cerr << "<<DEBUG>> End of gesvd" << std::endl;
        free(panel);
        dplasma_hqr_finalize( &qrtree );

        dague_ddesc_destroy((dague_ddesc_t*)matrixCI);
        dague_ddesc_destroy((dague_ddesc_t*)vectorCI);
        std::cerr << "<<DEBUG>> Data freed" << std::endl;

        free(matrixCI);
        free(vectorCI);
    }

    return rc;
}

void rational_fitter_parsec_one::fillci_p(int i0, int M, int N,
                                          double *CI, int ldci,
                                          const void *data,
                                          void* rptr)
{
    const vertical_segment* d = (const vertical_segment*)data;
    rational_function_1d *r = (rational_function_1d*)rptr;
    double *lCI;

    /* Let's start by dealing with the odd row */
    if ( i0%2 == 1 && M > 0 ) {
        lCI = CI;
        vec xi = d->get(i0 / 2);

        for(int j=0; j<N; ++j, lCI += ldci) {
            // Filling the p part
            const double pi = r->p(xi, j);
            lCI[0] = -pi;
        }
        M--;
        CI++;
        i0++;
    }

    int istart = i0 / 2;
    int iend   = (i0+M) / 2;

    for( int i=istart; i<iend; i++) {
        // Each constraint (fitting interval or point
        // add another dimension to the constraint
        // matrix

        lCI = CI + 2 * (i-istart);
        vec xi = d->get(i) ;

        // A row of the constraint matrix has this
        // form: [p_{0}(x_i), .., p_{np}(x_i)]
        // For the lower constraint and negated for
        // the upper constraint
        for(int j=0; j<N; ++j, lCI += ldci) {

            // Filling the p part
            const double pi = r->p(xi, j);

            lCI[0] =  pi;
            lCI[1] = -pi;
        }
    }

    /* Let's do last row if needed */
    if ( M%2 == 1  ) {
        lCI = CI + M - 1;
        vec xi = d->get(iend);

        for(int j=0; j<N; ++j, lCI += ldci) {
            // Filling the p part
            const double pi = r->p(xi, j);
            lCI[0] = -pi;
        }
    }
}

void rational_fitter_parsec_one::fillci_q(int i0, int M, int N,
                                          double *CI, int ldci,
                                          const void *data,
                                          void* rptr, int ny)
{
    const vertical_segment* d = (const vertical_segment*)data;
    rational_function_1d *r = (rational_function_1d*)rptr;
    double *lCI;

    /* Let's start by dealing with the odd row */
    if ( i0%2 == 1 && M > 0 ) {
        lCI = CI;
        vec xi = d->get(i0 / 2);

        for(int j=0; j<N; ++j, lCI += ldci) {
            vec yl, yu ;
            d->get(i0/2, yl, yu) ;

            const double qi = r->q(xi, j);

            lCI[0] = yl[ny] * qi;
        }
        M--;
        CI++;
        i0++;
    }

    int istart = i0 / 2;
    int iend   = (i0+M) / 2;

    for( int i=istart; i<iend; i++) {
        // Each constraint (fitting interval or point
        // add another dimension to the constraint
        // matrix

        lCI = CI + 2 * (i-istart);
        vec xi = d->get(i) ;

        // A row of the constraint matrix has this
        // form: [-f(x_i) q_{0}(x_i), .., -f(x_i) q_{nq}(x_i)]
        // For the lower constraint and negated for
        // the upper constraint
        for(int j=0; j<N; ++j, lCI += ldci) {

            vec yl, yu ;
            d->get(i, yl, yu) ;

            const double qi = r->q(xi, j);

            lCI[0] = -yu[ny] * qi;
            lCI[1] =  yl[ny] * qi;
        }
    }

    /* Let's do last row if needed */
    if ( M%2 == 1  ) {
        lCI = CI + M - 1;
        vec xi = d->get(iend);

        for(int j=0; j<N; ++j, lCI += ldci) {

            vec yl, yu ;
            d->get(iend, yl, yu) ;

            const double qi = r->q(xi, j);

            lCI[0] = -yu[ny] * qi;
        }
    }
}

int rational_fitter_parsec_one::solve_wrapper(int np, int nq, int M,
                                              double *ptr_CI, double *ptr_ci,
                                              void *rptr)
{
    rational_function_1d *r = (rational_function_1d*)rptr;
    const int N = np+nq ;

    // Compute the solution
    QuadProgPP::Matrix<double> CE(0.0, N, 0) ;
    QuadProgPP::Vector<double> ce(0.0, 0) ;
    QuadProgPP::Matrix<double> G (0.0, N, N) ;
    QuadProgPP::Vector<double> g (0.0, N) ;
    QuadProgPP::Vector<double> x;

    QuadProgPP::Matrix<double> CI(ptr_CI, N, M) ;
    QuadProgPP::Vector<double> ci(ptr_ci, M) ;

    // Select the size of the result vector to
    // be equal to the dimension of p + q
    for(int i=0; i<N; ++i) {
        G[i][i] = 1.0 ;
    }

    double cost = QuadProgPP::solve_quadprog(G, g, CE, ce, CI, ci, x);

    bool solves_qp = !(cost == std::numeric_limits<double>::infinity());
    for(int i=0; i<N; ++i) {
        const double v = x[i];
        solves_qp = solves_qp && !isnan(v)
            && (v != std::numeric_limits<double>::infinity()) ;
    }

    if(solves_qp) {
        // Recopy the vector d
        vec p(np), q(nq);
        double norm = 0.0 ;
        for(int i=0; i<N; ++i) {
            const double v = x[i];
            norm += v*v ;
            if(i < np) {
                p[i] = v ;
            }
            else {
                q[i - np] = v ;
            }
        }

        r->update(p, q);
#ifdef DEBUG
        std::cout << "<<INFO>> got solution " << *r << std::endl ;
#endif
        return norm > 0.0;
    }
    else {
        return 0;
    }
}
