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
            RUSIN_TH_TD,
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

        //! \brief static function for input type convertion. This
        //! function allocate the resulting vector.
        static double* convert(const double* invec, params::type intype,
                               params::type outtype)
        {
            int dim = dimension(outtype); // Get the size of the output vector

            if(dim > 0)
            {
                double* outvec = new double[dim];
                double  temvec[6]; // Temp CARTESIAN vectors
                to_cartesian(invec, intype, temvec);
                from_cartesian(temvec, outtype, outvec);

                return outvec;
            }
            else
            {
                return NULL;
            }
        }

        //! \brief static function for input type convertion. The outvec
        //! resulting vector should be allocated with the correct
        //! output size.
        static void convert(const double* invec, params::type intype,
                            params::type outtype, double* outvec)
        {
            double  temvec[6]; // Temp CARTESIAN vectors
            to_cartesian(invec, intype, temvec);
            from_cartesian(temvec, outtype, outvec);
        }

        //! \brief convert a input vector in a given parametrization to an
        //! output vector in a cartesian parametrization, that is two 3d
        //! vectors concatenated.
        static void to_cartesian(const double* invec, params::type intype,
                                 double* outvec)
        {
            double half[3];

            switch(intype)
            {
                // 2D Parametrizations
                case params::COS_TH_TD:
                    half_to_cartesian(acos(invec[0]), 0.0, acos(invec[1]), 0.0, outvec);
                    break;

                case params::RUSIN_TH_TD:
                    half_to_cartesian(invec[0], 0.0, invec[1], 0.0, outvec);
                    break;

                // 3D Parametrization
                case params::RUSIN_TH_PH_TD:
                    half_to_cartesian(invec[0], invec[1], invec[2], 0.0, outvec);
                    break;

                // 4D Parametrization
                case params::RUSIN_TH_PH_TD_PD:
                    half_to_cartesian(invec[0], invec[1], invec[2], invec[3], outvec);
                    break;


                // 6D Parametrization
                case params::CARTESIAN:
                    memcpy(outvec, invec, 6*sizeof(double));
                    break;

                default:
                    throw("Transformation not implemented, params::to_cartesian");
                    break;
            }

        }

        //! \brief convert a input CARTESIAN vector, that is two 3d vectors
        //! concatenated  to an output vector in a given parametrization.
        static void from_cartesian(const double* invec, params::type outtype,
                                   double* outvec)
        {
            // Compute the half vector
            double half[3] ;
            half[0] = invec[0] + invec[3];
            half[1] = invec[1] + invec[4];
            half[2] = invec[2] + invec[5];
            double half_norm = sqrt(half[0]*half[0] + half[1]*half[1] + half[2]*half[2]);
            half[0] /= half_norm;
            half[1] /= half_norm;
            half[2] /= half_norm;

            switch(outtype)
            {
                // 2D Parametrizations
                case params::COS_TH_TD:
                    outvec[0] = half[2];
                    outvec[1] = half[0]*outvec[0] + half[1]*outvec[1] + half[2]*outvec[2];
                    break;

                // 3D Parametrization
                case params::RUSIN_TH_PH_TD:
                    outvec[0] = acos(half[2]);
                    outvec[1] = atan2(half[0], half[1]);
                    outvec[2] = acos(half[0]*outvec[0] + half[1]*outvec[1] + half[2]*outvec[2]);
                    break;

                // 6D Parametrization
                case params::CARTESIAN:
                    memcpy(outvec, invec, 6*sizeof(double));
                    break;

                default:
                    throw("Transformation not implemented, params::from_cartesian");
                    break;
            }
        }

        //! \brief provide a dimension associated with a parametrization
        static int  dimension(params::type t)
        {
            switch(t)
            {
                // 2D Parametrizations
                case params::ISOTROPIC_TD_PD:
                case params::ROMEIRO_TH_TD:
                case params::COS_TH_TD:
                    return 2;
                    break;

                // 3D Parametrization
                case params::RUSIN_TH_PH_TD:
                case params::RUSIN_TH_TD_PD:
                case params::ISOTROPIC_TV_TL_DPHI:
                    return 3;
                    break;

                // 4D Parametrization
                case params::RUSIN_TH_PH_TD_PD:
                case params::SPHERICAL_TL_PL_TV_PV:
                    return 4;
                    break;

                // 6D Parametrization
                case params::CARTESIAN:
                    return 6;
                    break;

                default:
                    return -1;
                    break;
            }
        }

        //! \brief from the 4D definition of a half vector parametrization,
        //! export the cartesian coordinates.
        static void half_to_cartesian(double theta_h, double phi_h,
                                      double theta_d, double phi_d, double* out)
        {
            // Calculate the half vector
            double half[3];
            half[0] = sin(theta_h)*cos(phi_h);
            half[1] = sin(theta_h)*sin(phi_h);
            half[2] = cos(theta_h);

            // Compute the light vector using the rotation formula.
            out[0] = sin(theta_d)*cos(phi_d);
            out[1] = sin(theta_d)*sin(phi_d);
            out[2] = cos(theta_d);
            rotate_binormal(out, theta_h);
            rotate_normal(out, phi_h);

            // Compute the out vector from the in vector and the half
            // vector.
            const double dot = out[0]*half[0] + out[1]*half[1] + out[2]*half[2];
            out[3] = -out[0] + 2.0*dot * half[0];
            out[4] = -out[1] + 2.0*dot * half[1];
            out[5] = -out[2] + 2.0*dot * half[2];
        }

        //! \brief rotate a cartesian vector with respect to the normal of
        //! theta degrees.
        static void rotate_normal(double* vec, double theta)
        {
            const double cost = cos(theta);
            const double sint = sin(theta);

            vec[0] = cost * vec[0] - sint * vec[1];
            vec[1] = sint * vec[0] + cost * vec[1];
        }

        //! \brief rotate a cartesian vector with respect to the bi-normal of
        //! theta degrees.
        static void rotate_binormal(double* vec, double theta)
        {
            const double cost = cos(theta);
            const double sint = sin(theta);

            vec[0] = cost * vec[0] - sint * vec[2];
            vec[2] = sint * vec[0] + cost * vec[2];
        }
};
