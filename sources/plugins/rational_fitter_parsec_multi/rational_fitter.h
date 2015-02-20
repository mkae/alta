#pragma once

// Include STL
#include <vector>
#include <string>

// Interface
#include <core/function.h>
#include <core/rational_function.h>
#include <core/data.h>
#include <core/vertical_segment.h>
#include <core/fitter.h>
#include <core/args.h>
#include "gesvdm2.h"

/*! \brief A vertical segment fitter for rational function using the library QuadProg++
 *  You can find the library here: http://quadprog.sourceforge.net/
 *  \ingroup plugins
 */
class rational_fitter_parsec_multi : public fitter
{
  public: // methods

    rational_fitter_parsec_multi();
    virtual ~rational_fitter_parsec_multi();

    // Fitting a data object
    //
    virtual bool fit_data(const ptr<data>& d, ptr<function>& fit, const arguments& args);

    // Provide user parameters to the fitter
    //
    virtual void set_parameters(const arguments& args);

  private: // methods

    static void fill_p(const gesvdm2_args_t *args, int ny, int i0, int M,
		       double *P, int ldp);

    static void fill_q(const gesvdm2_args_t *args, int ny, int i0, int M,
		       double *P, int ldp);

    static void solve_init( const gesvdm2_args_t *pb, subproblem_t *spb );
    static int  solve_wrapper( const gesvdm2_args_t *args, subproblem_t *pb, int M, int ny,
			       double *CIptr, double *ciptr );
    static int solve_finalize( const gesvdm2_args_t *pb, subproblem_t *spb );

    static int  test_all_constraint( const vertical_segment     *data,
				     const rational_function_1d *r,
				     int ny );

    void initPbRf( subproblem_t *pb );

  protected: // data

    // Fitting a data object using np elements in the numerator and nq
    // elements in the denominator
    virtual bool fit_data(vertical_segment *d, int N, rational_function *rf, int &np);

    // min and Max usable np and nq values for the fitting
    int _max_np, _max_nq;
    int _min_np, _min_nq;

    // Add constraints to the boundary of the domain. You can shrink it of
    // the parameter --boundary-constraint *double*
    double _boundary;

    dague_context_t *_dague;
    int              _nbcores;
    const arguments *_args;

};
