/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2015 Université de Montréal

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#pragma once

// ALTA includes
#include "params.h"
#include "data.h"
#include "function.h"

// STL includes
#include <map>

/* \brief Error metrics class
 *
 * \details
 * This class contains the functionnalities to compute various error metric
 * between a function object and a data object or an interpolation data object
 * and a reference data object.
 *
 * Note: This class only contain static methods to regroup the different metric
 * and calling methods.
 *
 * \author Laurent Belcour <laurent.belcour@umontreal.ca>
 */
class errors {
   public:
      /* The result of an error metric is stored as a key/value combinaison here.
       * For example, the L2 error would be stored in {'L2': [float, ..., float]}.
       */
      typedef std::map<std::string, vec> metrics;

      /* Computes different metrics to compate two data objects 'in' to 'ref'.
      */
      static void compute(const data* in, const data* ref, metrics& res);

   private:

      /* Evaluate the output Eigen matrices containing the evaluation of the
       * interpolated data 'in' with respect to 'ref' abscissas.
       */
      static void evaluate(const data* in, const data* ref,
                           Eigen::MatrixXd& in_y, Eigen::MatrixXd& ref_y);

      /* Use Eigen Norm computation to compute efficiently a bunch of metrics
       * such a Lp norm, MSE and RMSE. The \a evaluate method needs to be used
       * to evaluate the 'in' and 'ref' arguments.
       */
      static void fastNormComputation(const Eigen::MatrixXd& in,
                                      const Eigen::MatrixXd& ref,
                                      Eigen::VectorXd& L1,
                                      Eigen::VectorXd& L2,
                                      Eigen::VectorXd& L3,
                                      Eigen::VectorXd& LInf,
                                      Eigen::VectorXd& mse,
                                      Eigen::VectorXd& rmse);
};
