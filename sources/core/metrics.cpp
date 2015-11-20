#include "metrics.h"

void errors::compute(const data* in, const data* ref, metrics& res) {
   //Here we go new way and faster because we evaluate the function just once
   Eigen::MatrixXd ref_y = Eigen::MatrixXd::Zero(ref->size(), ref->dimY());
   Eigen::MatrixXd inp_y = Eigen::MatrixXd::Zero(ref->size(), ref->dimY()) ;

#ifdef DEBUG
   timer  t;
   t.start();
#endif
   evaluate(in, ref, inp_y, ref_y);
#ifdef DEBUG
   t.stop();
   std::cout << "<<INFO>> Evaluate function for all data in  " << t << std::endl;
   t.reset();
#endif

   // Ouput norms
   Eigen::VectorXd L1_norm   = Eigen::VectorXd::Zero(ref->dimY());
   Eigen::VectorXd L2_norm   = Eigen::VectorXd::Zero(ref->dimY());
   Eigen::VectorXd L3_norm   = Eigen::VectorXd::Zero(ref->dimY());
   Eigen::VectorXd LInf_norm = Eigen::VectorXd::Zero(ref->dimY());
   Eigen::VectorXd mse       = Eigen::VectorXd::Zero(ref->dimY());
   Eigen::VectorXd rmse      = Eigen::VectorXd::Zero(ref->dimY());

#ifdef DEBUG
   t.start();
#endif
   fastNormComputation(inp_y, ref_y,
                       L1_norm, L2_norm, L3_norm,
                       LInf_norm, mse, rmse);
#ifdef DEBUG
   t.stop();
   std::cout << "<<INFO>> Fast Norm Computations in  " << t << std::endl;
#endif

   res.clear();
   res["L1"]   = L1_norm;
   res["L2"]   = L2_norm;
   res["L3"]   = L3_norm;
   res["LInf"] = LInf_norm;
   res["MSE"]  = mse;
   res["RMSE"] = rmse;
}


void errors::evaluate(const data* in,
                      const data* ref,
                      Eigen::MatrixXd& inp_y,
                      Eigen::MatrixXd& ref_y) {

   // TODO handle the case of output format conversion
   if(ref->output_parametrization() != in->output_parametrization()) {
      NOT_IMPLEMENTED();
   }

   // Temp variables
   vec ref_xy = vec::Zero(ref->dimX() + ref->dimY());
   vec ref_x  = vec::Zero(ref->dimX());
   vec dat_x  = vec::Zero(in->dimX());
   vec cart   = vec::Zero(6);

   // Constants
   const auto nY = ref->dimY();
   const auto nX = ref->dimX();

   // Evaluate the input data at each position of data_x configuration
   for(auto i=0; i<ref->size(); i++) {

      ref_xy = ref->get(i);
      ref_x  = ref_xy.head(nX);

      params::convert(ref_x.data(),
                      ref->input_parametrization(),
                      params::CARTESIAN,
                      cart.data());

      // Check if the output configuration is below the hemisphere when
      // converted to cartesian coordinates. Note that this prevent from
      // converting BTDF data.
      if(cart[2] >= 0.0 || cart[5] >= 0.0) {

         // Convert to the query data `dat` parametrization
         params::convert(cart.data(),
                         params::CARTESIAN,
                         in->input_parametrization(),
                         dat_x.data());

         ref_y.row(i) = ref_xy.tail(nY);
         params::convert(in->value(dat_x).data(),
                         in->output_parametrization(),
                         in->dimY(),
                         ref->output_parametrization(),
                         ref->dimY(),
                         inp_y.row(i).data());
      } else {
         ref_y.row(i).setZero();
         inp_y.row(i).setZero();
      }
   }
}

void errors::fastNormComputation(const Eigen::MatrixXd& o_data_y,
                                 const Eigen::MatrixXd& f_y,
                                 Eigen::VectorXd& L1,
                                 Eigen::VectorXd& L2,
                                 Eigen::VectorXd& L3,
                                 Eigen::VectorXd& LInf,
                                 Eigen::VectorXd& mse,
                                 Eigen::VectorXd& rmse) {
   assert(o_data_y.rows() == f_y.rows());
   assert(o_data_y.cols() == f_y.cols());
   const auto dMatrix = (o_data_y - f_y);

   // Check that the output dimensions match in size for all the
   // vectors and the difference matrix.
   const auto ncols = dMatrix.cols();
   assert(L1.size()   == ncols);
   assert(L2.size()   == ncols);
   assert(L3.size()   == ncols);
   assert(LInf.size() == ncols);

   // Compute the different distance metrics per output dimension
   for(auto i=0; i<ncols; i++) {
      L1(i)   = dMatrix.col(i).lpNorm<1>();
      L2(i)   = dMatrix.col(i).lpNorm<2>();
      L3(i)   = dMatrix.col(i).lpNorm<3>();
      LInf(i) = dMatrix.col(i).lpNorm<Eigen::Infinity>();
   }

   // Compute the RMSE and MSE
   mse = dMatrix.array().square().colwise().sum() / o_data_y.rows();
   rmse = mse.cwiseSqrt();
}
