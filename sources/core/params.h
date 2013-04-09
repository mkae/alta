#pragma once

/*! \brief a static class allowing to change from one parametrization
 *  to another.
 *  \ingroup core
 *
 *  \details
 *  Any \typedef function object or \typedef data object should have
 *  an associated parametrization.
 */
class params
{
    public: // data

        //! \brief list of all supported parametrization. An unsupported
        //! parametrization will go under the name <em>unknown</em>.
        enum type
        {
            ROMEIRO_TH_TD,
            RUSIN_TH_PH_TD,
            RUSIN_TH_TD_PD,
            RUSIN_TH_PH_TD_PD,
            COS_TH_TD,
            ISOTROPIC_TV_TL_DPHI,
            ISOTROPIC_TD_PD, // Difference between two directions such as R and H
            CARTESIAN,
            SPHERICAL_TL_PL_TV_PV,
            UNKNOWN
        };

    public: // methods

        //! \brief static function for input type convertion.
        //! \todo implement the convertion method using a case on the
        //! out type. It should work if the input vector is converted
        //! to a 6D vector (in, out).
        static vec convert(const vec& invec, params::type intype, params::type outtype)
        {
            return invec;
        }
};
