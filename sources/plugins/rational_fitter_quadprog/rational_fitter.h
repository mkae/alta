/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014 Inria

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
 *  \ingroup fitters
 *  \todo : WRITE MORE DOCUMENTATION
 */
class rational_fitter_quadprog : public fitter
{
  public: // methods

    rational_fitter_quadprog() ;
    virtual ~rational_fitter_quadprog() ;

    // Fitting a data object
    //
    virtual bool fit_data(const ptr<data>& d, ptr<function>& fit, const arguments& args) ;

    // Provide user parameters to the fitter
    //
    virtual void set_parameters(const arguments& args) ;

  private: // methods

  protected: // data

    // Fitting a data object using np elements in the numerator and nq
    // elements in the denominator
    virtual bool fit_data(const ptr<vertical_segment>& d, int np, int nq, const ptr<rational_function>& fit) ;
    virtual bool fit_data(const ptr<vertical_segment>& dat, int np, int nq, int ny, rational_function_1d* fit) ;

    // min and Max usable np and nq values for the fitting
    int _max_np, _max_nq ;
    int _min_np, _min_nq ;

    // Add constraints to the boundary of the domain. You can shrink it of
    // the parameter --boundary-constraint *double*
    double _boundary;

    std::string _scheduling_mode;
    int _scheduling_chunk_size;
    double _scheduling_grow_factor;
    bool _export_qp;
    double _delta;
    bool _add_ls_energy;
} ;

