/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2014 Bordeaux-INP
   Copyright (C) 2014 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

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

/*! \brief A vertical segment fitter for rational function using the library QuadProg++
 *  You can find the library here: http://quadprog.sourceforge.net/
 *  \ingroup plugins
 */
class rational_fitter_parsec_one : public fitter
{
  private:
    dague_context_t *dague;
    int              nbcores;

  public: // methods

    rational_fitter_parsec_one() ;
    virtual ~rational_fitter_parsec_one() ;

    // Fitting a data object
    //
    virtual bool fit_data(const ptr<data> d, function* fit, const arguments& args) ;

    // Provide user parameters to the fitter
    //
    virtual void set_parameters(const arguments& args) ;

    static void fillci_p(int i0, int M, int N,
                         double *CI, int ldci,
                         const void *data,
                         void* rptr);
    static void fillci_q(int i0, int M, int N,
                         double *CI, int ldci,
                         const void *data,
                         void* rptr, int ny);

    static int solve_wrapper(int np, int nq, int M,
                             double *ptr_CI, double *ptr_ci,
                             void *rptr);

  private: // methods

  protected: // data

    // Fitting a data object using np elements in the numerator and nq
    // elements in the denominator
    virtual bool fit_data(const vertical_segment* d, int np, int nq, rational_function* fit) ;
    virtual bool fit_data(const vertical_segment* dat, int np, int nq, int ny, rational_function_1d* fit) ;

    // min and Max usable np and nq values for the fitting
    int _max_np, _max_nq ;
    int _min_np, _min_nq ;

    // Add constraints to the boundary of the domain. You can shrink it of
    // the parameter --boundary-constraint *double*
    double _boundary;
} ;
