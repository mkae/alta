#pragma once

#include <Eigen/SVD>
#include <Array.hh>
#include <QuadProg++.hh>

#include <core/rational_function.h>
#include <core/vertical_segment.h>

class quadratic_program
{
public:
    //! \brief Constructor need to specify the number of coefficients
    quadratic_program(int np, int nq) : _np(np), _nq(nq), CI(0.0, _np+_nq, 0) { }

    //! \brief Remove the already defined constraints
    void clear_constraints()
    {
        CI.resize(_np+_nq, 0);
    }

    //! \brief Add a constraint by specifying the vector
    void add_constraints(const QuadProgPP::Vector<double> vec)
    {
        const int m = CI.nrows();
        const int n = CI.ncols();

        if(n > 0)
        {
            // Construct temp buffer
            double* temp = new double[n*m];
            for(int u=0; u<n; ++u)
            {
                for(int v=0; v<m; ++v)
                {
                    temp[u*m + v] = CI[v][u];
                }
            }

            // Resize matrix CI
            CI.resize(m, n+1);

            // Recopy data
            for(int u=0; u<n+1; ++u)
            {
                if(u==n)
                {
                    for(int v=0; v<m; ++v)
                        CI[v][u] = vec[v];
                }
                else
                {
                    for(int v=0; v<m; ++v)
                        CI[v][u] = temp[u*m + v];
                }
            }
            delete[] temp;
        }
        else
        {
            // Resize matrix CI
            CI.resize(m, 1);

            // Recopy data
            for(int u=0; u<m; ++u)
                CI[n][u] = vec[u];
        }
    }

    //! \brief Solves the quadratic program
    bool solve_program(QuadProgPP::Vector<double>& v, double& delta)
    {        
        const int m = CI.nrows();
        const int n = CI.ncols();

        QuadProgPP::Matrix<double> G (0.0, m, m) ;
        QuadProgPP::Vector<double> g (0.0, m) ;
        QuadProgPP::Vector<double> ci(0.0, n) ;
        QuadProgPP::Matrix<double> CE(0.0, m, 0) ;
        QuadProgPP::Vector<double> ce(0.0, 0) ;

        // Update the ci column with the delta parameter
        // (See Celis et al. 2007 p.12)
        Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::HouseholderQRPreconditioner> svd(Eigen::MatrixXd::Map(&CI[0][0], m, n));
        const double sigma_m = svd.singularValues()(std::min(m, n)-1) ;
        const double sigma_M = svd.singularValues()(0) ;
        delta = sigma_M / sigma_m ;

        // Select the size of the result vector to
        // be equal to the dimension of p + q
        for(int i=0; i<m; ++i)
        {
            G[i][i] = 1.0 ;
        }

        // Each constraint (fitting interval or point
        // add another dimension to the constraint
        // matrix
        for(int i=0; i<n; ++i)
        {
            // Norm of the row vector
            double norm = 0.0 ;

            for(int j=0; j<m; ++j)
            {
                norm += CI[j][i]*CI[j][i] ; ;
            }

            // Set the c vector, will later be updated using the
            // delta parameter.
            ci[i] = -delta * sqrt(norm) ;
        }

        // Compute the solution
        const double cost = QuadProgPP::solve_quadprog(G, g, CE, ce, CI, ci, v);

        bool solves_qp = !(cost == std::numeric_limits<double>::infinity());
        return solves_qp;
    }

    //! \brief Test all the constraints of the data
    bool test_constraints(const rational_function* r, const vertical_segment* data)
    {
        int nb_failed = 0;
        for(int n=0; n<data->size(); ++n)
        {
            vec x, yl, yu;
            data->get(n, x, yl, yu);

            vec y = r->value(x);
            if(y < yl || y > yu)
            {
                nb_failed++;
            }
        }

#ifdef DEBUG
        std::cout << "<<TRACE>> " << nb_failed << " constraints where not satified." << std::endl;
#endif

        return nb_failed == 0;
    }

protected:
    int _np, _nq ;
    QuadProgPP::Matrix<double> CI;

};
