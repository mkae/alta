/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014, 2015, 2016 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#pragma once

// Include STL
#include <vector>
#include <string>
#include <iostream>

// Interface
#include "common.h"
#include "function.h"
#include "data.h"
#include "fitter.h"
#include "args.h"

namespace alta {

/*! \ingroup core
 *  \ingroup datas
 *
 *  \brief
 *  A vertical segment data class
 *
 *  This class implement a data representation of vertical segments in the
 *  sense of Pacanowski et al. [2012]. Each data point is in fact composed
 *  of a middle point \f$ x \f$ and an upper \f$ \overline{x} \f$ and lower
 *  bound \f$ \underline{x} \f$.
 *
 *  To retrieve the complete vertical segment data \f$ [x, \underline{x},
 *  \overline{x}] \f$, a special function is provided. The functions
 *  inherited from \a data will only return the middle point.
 *
 *  It is possible to load regular ALTA file using a vertical segment data
 *  loader. It will automatically generate vertical segments. You can
 *  control the behavior of the vertical segments using the following
 *  option in the command line:
 *
 *  + <strong>\-\-dt</strong> specify the size of the vertical segment. If the
 *     option <strong>\-\-dt-relative</strong> is not set, this size is absolute: \f$ [x,
 *     x - dt, x + dt] \f$. If the <strong>\-\-dt-relative</strong> option is set, the
 *     vertical segment size is relative to the middle point value \f$ x \f$:
 *     \f$ [x, x (1 - dt), x (1 + dt)] \f$. You can specify the vertical
 *     segment to be equal to the max of the relative and absolute sizes
 *     using the <strong>\-\-dt-max</strong> option.
 *
 *   + <strong>\-\-data-positive</strong> for the vertical segment to stay in the
 *    positive region. The negative values are replaced by zeros.
 *
 *
 *  The data of the vertical segment can be restricted to subpart of the
 *  original data by specifying the bounding box of the input and output
 *  domain:
 *
 *  + <strong>\-\-min</strong> <em>[vec]</em> specify the minimun input
 *    coordinate that should be loaded. All data with input coordinate
 *    less than this vector will be discarded.
 *
 *  + <strong>\-\-max</strong> <em>[vec]</em> specify the maximum input
 *    coordinate that should be loaded. All data with input coordinate
 *     greater than this vector will be discarded.
 *
 *   + <strong>\-\-ymin</strong> <em>[vec]</em> specify the minimun output
 *    coordinate that should be loaded. All data with associated value
 *    less than this vector will be discarded.
 *
 *   + <strong>\-\-ymax</strong> <em>[vec]</em> specify the maximum output
 *    coordinate that should be loaded. All data with associated value
 *    greater than this vector will be discarded.
 */
class vertical_segment : public data
{
   public: // methods

      vertical_segment(const parameters& params,
                       size_t size,
                       std::shared_ptr<double> data);

      vertical_segment(const parameters& params, unsigned int size)
          ALTA_DEPRECATED;

      // Acces to data
      virtual vec get(int i) const ;

      virtual vec value(const vec&) const {
         NOT_IMPLEMENTED();
      }

      //! \brief Put the sample inside the data at index I.
      virtual void set(int i, const vec& x);

      //! \brief Specific accessor to a vertical segment, this gives the
      //! complete vector, plus the ordinate segment
      virtual void get(int i, vec &x, vec &yl, vec &yu) const ;

      //! \brief Specific accessor to a vertical segment. Provides only the
      //! ordinate segment.
      virtual void get(int i, vec& yl, vec& yu) const ;

   private: // method

      //! \brief Return the number of columns in each row.
      size_t column_number() const
      {
          return _parameters.dimX() + 3 * _parameters.dimY();
      }

      // Return a matrix view of the data: all the Xi and Yi followed by
      // confidence interval data (lower and upper bound of the Yi).  Thus,
      // it has (dimX + 3 * dimY) columns and SIZE rows.
      Eigen::Map<Eigen::MatrixXd> matrix_view() const
      {
          return Eigen::Map<Eigen::MatrixXd>(_data.get(), size(),
                                             column_number());
      }

      // Return a matrix view of DATA that excludes confidence interval data.
      // It has (dimX + dimY) columns and SIZE rows.
      Eigen::Map<Eigen::MatrixXd, 0, Eigen::OuterStride<> > data_view() const
      {
          return Eigen::Map<Eigen::MatrixXd, 0, Eigen::OuterStride<> >
              (_data.get(), _parameters.dimX() + _parameters.dimY(), size(),
               Eigen::OuterStride<>(column_number()));
      }

  protected: // method

      //! \brief From a correct input configuration 'x' with size
      //! dimX()+dimY(), generate a vertical ! segment satisfying this object's
      //! parameters.
      virtual vec vs(const vec& x) const;

	protected: // data

      // Store for each point of data, the upper and lower value.
      std::shared_ptr<double> _data;

      // Store the different arguments for the vertical segment: is it using
      // relative or absolute intervals? What is the dt used ?
      bool   _is_absolute;
      double _dt;
};
}

/* -*- c++ -*- */
