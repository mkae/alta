/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2015 Inria

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

// ALTA includes
#include <core/params.h>
#include <tests.h>

using namespace alta::tests;

#include <cppqc.h>                                // CppQuickCheck

using namespace cppqc;

// STL includes
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cstdlib>


typedef double Angle;

// Trivial adaptation of the 'ChooseStatelessGenerator' template class that
// targets integers.
template<class Real>
class ChooseStatelessGenerator
{
public:
    ChooseStatelessGenerator(Real min, Real max): m_min(min), m_max(max)
    {
        assert(min <= max);
    }

    Real unGen(RngEngine &rng, std::size_t) const
    {
        // XXX: We'd like to get a non-uniform distribution, with more
        // entries close to M_MIN and M_MAX.
        boost::uniform_real<Real> dist(m_min, m_max);
        return dist(rng);
    }

    std::vector<Real> shrink(Real x) const
    {
        std::vector<Real> ret;
        ret.reserve(x - m_min);
        if (abs(m_min) <= abs(m_max)) {
            for (Real i = m_min; i != x; ++i)
                ret.push_back(i);
        } else {
            for (Real i = m_max; i != x; --i)
                ret.push_back(i);
        }
        return ret;
    }

private:
    const Real m_min;
    const Real m_max;
};

/// Generates a random real in the range min..max inclusive, requires that
/// min <= max. Shrinks towards smaller absolute values.
template<class Real>
StatelessGenerator<Real > chooseReal(Real min, Real max)
{
    return ChooseStatelessGenerator<Real>(min, max);
}


// Cartesian-to-Rusinkiewicz valid conversion property with ωi co-linear with
// ωo.

// Whether to dump details of each test.
static bool verbose;

template<typename Angle>
struct PropCartToRusin: Property<Angle, Angle>
{
    PropCartToRusin():
        Property<Angle, Angle>(chooseReal(0., M_PI_2),
                               chooseReal(0., M_PI))
        {}

    bool check(const Angle& theta, const Angle& phi) const
    {
        auto v = spherical_to_cartesian(theta, phi);
        vec cart(6);
        cart[0] = v[0], cart[1] = v[1], cart[2] = v[2];
        cart[3] = v[0], cart[4] = v[1], cart[5] = v[2];

        vec x(4);
        params::convert(&cart[0], params::CARTESIAN,
                        params::RUSIN_TH_PH_TD_PD, &x[0]);

        auto theta_h = x[0];
        auto phi_h = x[1];
        auto theta_d = x[2];
        auto phi_d = x[3];

        // Make sure θh = θi and φh = φi, and θd = φd = 0.
        // Note that we have less precision on φd.
        bool pass = close_to(theta_h, (double)theta)
            && close_to(phi_h, (double)phi)
            && close_to(theta_d, 0.)
            && close_to(phi_d, 0., 5e-6);

        if (verbose || !pass)
            std::cout << std::setprecision(20)
                      << "conf " << theta << " " << phi
                      << " aka. cart = " << cart
                      << " → θh = " << theta_h
                      << " → φh = " << phi_h
                      << " → θd = " << theta_d
                      << " → φd = " << phi_d
                      << std::endl;

        return pass;
    }
    std::string name() const
    {
        return "Cartesian → 4D Rusinkiewicz conversion works";
    }
    std::string classify(const Angle& theta, const Angle& phi) const
    {
        if (close_to(theta, 0., 1e-2))
            return "grazing angles";
        else if (close_to(theta, M_PI_2, 1e-2))
            return "polar angles";
        else
            return "random angles";
    }
    bool trivial(const Angle& theta, const Angle& phi) const
    {
        return close_to(theta, 0.) && close_to(phi, 0.);
    }
};

int main(int argc, char** argv)
{
    verbose = getenv("VERBOSE") != NULL;
    auto status = quickCheckOutput(PropCartToRusin<Angle>(), std::cout, 200);
    return status.result == QC_SUCCESS ? EXIT_SUCCESS : EXIT_FAILURE;
}
