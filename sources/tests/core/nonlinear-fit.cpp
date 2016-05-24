/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2016 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

#include <core/data.h>
#include <core/vertical_segment.h>
#include <core/function.h>
#include <core/fitter.h>
#include <core/plugins_manager.h>
#include <tests.h>

#include <cstdlib>
#include <iostream>
#include <string>
#include <list>

using namespace alta;
using namespace alta::tests;

// Load input data, then iterate over nonlinear fitters and fit it.
int main(int argc, char *argv[])
{
    static std::list<std::string> fitters
        { "nonlinear_fitter_eigen",

          // The following fitters are optional.
          "nonlinear_fitter_nlopt",
          "nonlinear_fitter_ipopt",
          "nonlinear_fitter_ceres"
        };

    std::string input_file;
    {
        // Process the default test file.
        static const std::string data_file = "Kirby2.dat";
        std::string data_dir = getenv("TEST_DATA_DIRECTORY") != NULL
            ? getenv("TEST_DATA_DIRECTORY") : ".";
        input_file = data_dir + "/" + data_file;
    }

    auto data = ptr<alta::data>(new vertical_segment());
    data->load(input_file);

    for (auto&& fitter_name: fitters) {
        auto function = plugins_manager::get_function("nonlinear_function_diffuse",
                                                      data->parametrization());
        TEST_ASSERT(function != NULL);

        auto fitter = plugins_manager::get_fitter(fitter_name);

        // If FITTER_NAME denotes an optional fitter, FITTER may be NULL.
        // Just skip it.
        if (fitter != NULL) {
            std::cerr << "testing fitter '" << fitter_name << "'...\n";
            TEST_ASSERT(fitter->fit_data(data, function, arguments()));

            // Verify basic properties of FUNCTION.
            // XXX: "nonlinear_function_diffuse" always has dimX = 6.
            TEST_ASSERT(function->parametrization().dimX()
                        >= data->parametrization().dimX());
            TEST_ASSERT(function->parametrization().dimY()
                        >= data->parametrization().dimY());
        } else {
            std::cerr << "skipping fitter '" << fitter_name << "'\n";
        }
    }

    return EXIT_SUCCESS;
}
