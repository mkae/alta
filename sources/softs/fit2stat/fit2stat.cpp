/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2015 CNRS

   This file is part of ALTA.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0.  If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.  */

/*! \package fit2stats
 *  \ingroup commands
 *  \brief
 *  This command takes at input the result of a BRDF fit as well as data and outputs statistics regarding
 *  the fitting. 
 *
 *  \details
 *  
 *  <h3>Parameters</h3>
 *  <ul>
 */

//STL
#include <cstdlib>


//ALTA
#include <core/common.h>
#include <core/args.h>
#include <core/data.h>
#include <core/vertical_segment.h>
#include <core/function.h>
#include <core/fitter.h>
#include <core/plugins_manager.h>

using namespace alta;

//TODO: Move to ALTA CORE
//TODO: Documentation
//TODO:  Add the relative error computation
class Error
{
public:
  static vec meanSquareError(ptr<data> const & data, ptr<function> const & f )
  {
    vec mse = vec::Zero( f->dimY() );

    vec dat = vec::Zero( data->parametrization().dimX() + data->parametrization().dimY() );
    
    //Note that: data_x = dat.head( data->parametrization().dimX() );
    //           data_y = dat.tail( data->parametrization().dimY() );
    vec f_y = vec::Zero( f->dimY() );
    
    for(unsigned int i=0; i < (unsigned int) data->size(); i++)
    {
      dat = data->get(i); 

      f_y = f->value(dat.head( data->parametrization().dimX() ) );

      mse += (f_y - dat.tail( data->parametrization().dimY() ) ).array().square().matrix();
    }
    
    return mse * (1.0 / data->size() );

  }

  static vec rootMeanSquareError(ptr<data> const & data, ptr<function> const & f )
  {
    vec const rmse = meanSquareError(data, f);
    return rmse.cwiseSqrt();
  }

};

 
// TODO : Write the documentation
// TODO: Move this class to ALTA Core
// The methods provided by this class assume that the parametrization of 
// the data and function are the same.
class Norm
{
  public:
  static vec L2( ptr<data> const & data, ptr<function> const & f )
  {
    vec l2 = vec::Zero( f->dimY() );

    vec dat = vec::Zero( data->parametrization().dimX() + data->parametrization().dimY() );
    
    //Note that: data_x = dat.head( data->parametrization().dimX() );
    //           data_y = dat.tail( data->parametrization().dimY() );
    vec f_y = vec::Zero( f->dimY() );
    
    for(unsigned int i=0; i < (unsigned int) data->size(); i++)
    {
      dat = data->get(i);

      f_y = f->value(dat.head( data->parametrization().dimX() ) );
      
      l2 += ( f_y - dat.tail( data->parametrization().dimY() )).array().square().matrix();

    }
    return l2.cwiseSqrt();
  }

  static vec weightedL2( ptr<data> const & data, ptr<function> const & f, vec const & weights)
  {
    vec l2 = vec::Zero( f->dimY() );

    vec dat = vec::Zero( data->parametrization().dimX() + data->parametrization().dimY() );

    for(unsigned int i=0; i < (unsigned int) data->size(); i++)
    {
      dat = data->get(i);
      l2 += ((f->value(dat.head( data->parametrization().dimX() )) - dat.tail( data->parametrization().dimY() )) * weights[i]).array().square().matrix();

    }
    return l2.cwiseSqrt();
  }

  static vec LInf( ptr<data> const & data, ptr<function> const & f )
  {
    vec max = vec::Zero( f->dimY() );

    vec dat = vec::Zero( data->parametrization().dimX() + data->parametrization().dimY() );
  
    for(unsigned int i=0; i < data->size(); i++)
    {
      dat = data->get(i);          
      //vec const data_x = dat.head( data->parametrization().dimX() );
      //vec const data_y = dat.tail( data->parametrization().dimY() );

      max =  max.cwiseMax( (f->value(dat.head( data->parametrization().dimX() )) - dat.tail( data->parametrization().dimY() )).cwiseAbs() );

    }
    return max;
  }

  static vec L1( ptr<data> const & data, ptr<function> const & f )
  {
    vec L1 = vec::Zero( f->dimY() );

    vec dat = vec::Zero( data->parametrization().dimX() + data->parametrization().dimY() );

    for(unsigned int i=0; i < (unsigned int ) data->size(); i++)
    {
      dat = data->get(i);          
      //vec const data_x = dat.head( data->parametrization().dimX() );
      //vec const data_y = dat.tail( data->parametrization().dimY() );

      L1 += (f->value(dat.head( data->parametrization().dimX() )) - dat.tail( data->parametrization().dimY() )).cwiseAbs();

    }
    return L1;
  }

  /*! 
   * \brief Compute the Lp norm from a given p (with \f$ p >= 1 \f$) 
   * \f$ \big(\sum_{i=1}^{i=n} |x_i|^p\big)^{1/p} \f$
   */
  static vec Lp(ptr<data> const & data, ptr<function> const & f, double p)
  {
    //if p is not superior or equal to 1.0 this formula does not define a norm
    assert( p >= 1.0 );

    // This is another strategy to compute the norm
    // First we put all the data into Eigen::Array
    // Then we perform the computation on the big array directly
    Eigen::ArrayXd  a_Lp = Eigen::ArrayXd::Zero( f->dimY() );    
    Eigen::ArrayXXd  all_data_y( data->size(), f->dimY() );
    all_data_y.setZero(data->size(), f->dimY());

    Eigen::ArrayXXd  all_y_fx( data->size(), f->dimY() );
    all_y_fx.setZero(data->size(), f->dimY());

    //timer chrono;
    //chrono.start();
    for(unsigned int i=0; i < (unsigned int) data->size(); i++)
    {
      Eigen::ArrayXd const dat = data->get(i);          
      Eigen::ArrayXd const data_x = dat.head( data->parametrization().dimX() );
      Eigen::ArrayXd const data_y = dat.tail( data->parametrization().dimY() );
      Eigen::ArrayXd const y_fx = f->value(data_x);
      
      all_data_y.row(i) = data_y;
      all_y_fx.row(i)   = y_fx;

      //a_Lp += pow( abs(y_fx - data_y), p );
    }
    // chrono.stop();
    // std::cout << "<<INFO>> Time to demultiplex data " << chrono << std::endl;
    // chrono.reset();

    Eigen::ArrayXd pre_cal = pow(abs(all_data_y - all_y_fx), p).colwise().sum();


    double const one_over_p = 1.0 / p;

    return pre_cal.pow( one_over_p );

    //return Lp.array().pow( one_over_p );
    
    //Implici conversion from Eigen::Array to Matrix
    return a_Lp.pow( one_over_p );

  }

};

//Returns true if the conversion was successful
bool convertDataToFunctionParam(ptr<data> const & data,  
                                ptr<function> const & f,
                                bool & conversion_necessary,
                                vertical_segment* & converted_data )
{
  conversion_necessary = true;
  if( data->parametrization().input_parametrization() == f->input_parametrization() )
  {
    if( data->parametrization().output_parametrization() == f->output_parametrization() )
    {
      conversion_necessary = false;
      return true;
    }
    else // Ouput parametrizations are different. 
         // Output parametrization of  the Function prevails
    {
      converted_data = new vertical_segment( data->parametrization().input_parametrization(), 
                                             f->output_parametrization(),
                                             data->size() );
      parameters p(f->input_parametrization(),
                   f->output_parametrization());
      p.setDimX(f->dimX());
      p.setDimY(f->dimY());
      converted_data->setParametrization(p);

      //Note that: data_x = dat.head( data->parametrization().dimX() );
      //           data_y = dat.tail( data->parametrization().dimY() );
      vec dat        = vec::Zero( data->parametrization().dimX() + data->parametrization().dimY() );
      vec data_x     = vec::Zero( data->parametrization().dimX() );
      vec data_y     = vec::Zero( data->parametrization().dimY() );
      vec new_data_y = vec::Zero( f->dimY() );
      vec new_data   = vec::Zero(data->parametrization().dimX() + f->dimY() );

      for( unsigned int i=0; i < converted_data->size(); i++)
      {
        dat = data->get(i);          
        data_x = dat.head( data->parametrization().dimX() );
        data_y = dat.tail( data->parametrization().dimY() );

        params::convert(&data_y[0],
                        data->parametrization().output_parametrization(),
                        data->parametrization().dimY(),
                        f->output_parametrization(), f->dimY(),  &new_data_y[0]);

        new_data.head( data->parametrization().dimX() ) = data_x;
        new_data.tail( f->dimY() )    = new_data_y;

        converted_data->set(i, new_data );
      }
      
    }
  }// End if 
  else // Input parametrization are different
  {
    // Output param are the same
    // Converting to the function input parametrization
    if( data->parametrization().output_parametrization() == f->output_parametrization() ) 
    {
      parameters p(f->dimX(), f->dimY());
      p.setParametrization(f->input_parametrization());
      p.setParametrization(f->output_parametrization());

      converted_data = new vertical_segment( f->dimX(), f->dimY(), data->size() );
      converted_data->setParametrization(p);
      
      //Note that: data_x = dat.head( data->parametrization().dimX() );
      //           data_y = dat.tail( data->parametrization().dimY() );
      vec dat        = vec::Zero( data->parametrization().dimX() + data->parametrization().dimY() );
      vec data_x     = vec::Zero( data->parametrization().dimX() );
      vec data_y     = vec::Zero( data->parametrization().dimY() );
      
      vec new_data = vec::Zero( f->dimX() + f->dimY() );
      vec new_data_y = vec::Zero( f->dimY() );
      vec new_data_x = vec::Zero( f->dimX() );

      for( unsigned int i=0; i < (unsigned int) converted_data->size(); i++)
      {

        dat = data->get(i);
        data_x = dat.head( data->parametrization().dimX() );
        data_y = dat.tail( data->parametrization().dimY() );
        
        params::convert(&data_x[0],
                        data->parametrization().input_parametrization(),
                        f->input_parametrization(),
                        &new_data_x[0]);

        new_data.head( f->dimX() ) = new_data_x;
        new_data.tail( f->dimY() ) = data_y;

        converted_data->set(i, new_data );
      }//end of for-loop

    }
    else //EVERYTHING IS DIFFERENT ! ALL TO FUNCTION PARAMETRIZATIONS
    {
      //       unsigned int const new_size = f->dimX() + f->dimY();
      // converted_data = new vertical_segment( f->dimX(), f->dimY(), data->size() );
      // converted_data->setParametrizations( f->input_parametrization(), f->output_parametrization() );
      
      // for( unsigned int i=0; i < converted_data->size(); i++)
      // {
      //   vec const dat = data->get(i);
      //   vec const data_x = dat.head( data->parametrization().dimX() );
      //   vec const data_y = dat.tail( data->parametrization().dimY() );
      // }
      NOT_IMPLEMENTED();
    }
  }

  return true;
}


// TODO : Add a new parametriztion to params.h for the COS_TL_TV
//
void computeCosineFactorsFromData( ptr<data> const & data, 
                                   vec & cosine_theta_light, 
                                   vec & cosine_theta_view )
{
  cosine_theta_view.setZero( data->size() );
  cosine_theta_light.setZero( data->size() );
  
  if( data->parametrization().input_parametrization() == params::CARTESIAN)
  {
    vec data_x = vec::Zero( data->parametrization().dimX() );

    for( unsigned int i=0; i < data->size(); i++)
    {
      data_x = (data->get(i)).head( data->parametrization().dimX() );
      cosine_theta_light[i] = data_x[2];
      cosine_theta_view[i]  = data_x[5];
    }
  }
  else
  {
    NOT_IMPLEMENTED();
  }
}


void  demultiplexData( ptr<data> const & data, 
                       Eigen::ArrayXXd & o_data_x, 
                       Eigen::ArrayXXd & o_data_y )
{

  vec dat        = vec::Zero( data->parametrization().dimX() + data->parametrization().dimY() );
  vec data_x     = vec::Zero( data->parametrization().dimX() );
  vec data_y     = vec::Zero( data->parametrization().dimY() );
      
  for( unsigned int i=0; i < data->size(); i++)
  {
    dat    = data->get(i);
    data_x = dat.head( data->parametrization().dimX() );
    data_y = dat.tail( data->parametrization().dimY() );

    o_data_x.row(i) = data_x;
    o_data_y.row(i) = data_y;
  }

}

void evaluateFunctionAtData( Eigen::ArrayXXd const & data_x, 
                             ptr<function> const & f,
                             Eigen::ArrayXXd & f_y )
{
//  Eigen::VectorXd vdata_x
  for( unsigned int i=0; i < data_x.rows(); i++)
  {
      f_y.row(i) = f->value( data_x.row(i) ) ;
  }
}


void fastNormComputation( Eigen::ArrayXXd const & o_data_y, 
                          Eigen::ArrayXXd const & f_y,
                          Eigen::VectorXd & L1,
                          Eigen::VectorXd & L2, 
                          Eigen::VectorXd & L3, 
                          Eigen::VectorXd & LInf,
                          Eigen::VectorXd & mse,
                          Eigen::VectorXd & rmse )
{
  Eigen::MatrixXd  const distance_matrix = (o_data_y - f_y).matrix();

  for( unsigned int i=0; i < f_y.cols(); i++)
  {
    L1(i)   = distance_matrix.col(i).lpNorm<1>();
    L2(i)   = distance_matrix.col(i).lpNorm<2>();
    L3(i)   = distance_matrix.col(i).lpNorm<3>();
    LInf(i) = distance_matrix.col(i).lpNorm<Eigen::Infinity>();
  }

  mse = distance_matrix.array().square().colwise().sum() / o_data_y.rows();
  rmse = mse.cwiseSqrt();

}

void fastWeightedErrors( Eigen::ArrayXXd const & o_data_y, 
                         Eigen::ArrayXXd const & f_y,
                         Eigen::VectorXd const & weights,
                         Eigen::VectorXd & mse,
                         Eigen::VectorXd & rmse )
{

  Eigen::MatrixXd  const distance_matrix = (o_data_y - f_y).matrix();
  Eigen::MatrixXd tmp = Eigen::MatrixXd( distance_matrix.rows(), distance_matrix.cols() );

  for( unsigned int i=0; i < f_y.cols() ; i++)
  {
    tmp.col(i) = distance_matrix.col(i).cwiseProduct( weights );
  }

  mse = tmp.array().square().colwise().sum();
  
  mse /= o_data_y.rows();
  rmse = mse.cwiseSqrt();

}





int 
main(int argc, char* argv[])
{

  arguments args(argc, argv) ;

  if(args.is_defined("help") || args.is_defined("h")) 
  {
    std::cout << "Usage: fit2stat [options] --input brdf.file data.file " << std::endl ;
    std::cout << std::endl;
    std::cout << "Mandatory arguments:" << std::endl;
    std::cout << "  --brdf    [brdf_filename] : needs to be an alta file describing a BRDF" << std::endl;
    std::cout << "  --input   [data_filename] : needs to be an alta file describing  data" << std::endl;
    std::cout << std::endl;

    std::cout << "Optional arguments:" << std::endl;
    std::cout << "  --ymin     " << std::endl ;
    std::cout << "  --xmin     " << std::endl ;
    std::cout << "  --xmax    " << std::endl ;

    return EXIT_SUCCESS;
  }

  std::cout << "<<INFO>> Total System Memory:" 
            << plugins_manager::get_system_memory()  / (1024*1024*1024.0) << " Gigabytes "
            << std::endl;


  if(! args.is_defined("brdf")) 
  {
    std::cerr << "<<ERROR>> the BRDF filename is not defined" << std::endl ;
    return EXIT_FAILURE ;
  }

  if(! args.is_defined("input")) 
  {
    std::cerr << "<<ERROR>> the data filename is not defined" << std::endl ;
    return EXIT_FAILURE ;
  }

  std::cout << "<<INFO>> Loading data ..." << std::endl;
  ptr<vertical_segment> vs_data = ptr<vertical_segment>(new vertical_segment());
  //ptr<data> vs_data = plugins_manager::get_data("vertical_segment");
  
  timer  t;
  try
  {
    t.start();
    vs_data->load(args["input"], args);  
      t.stop();
  }
  catch(...)
  {
    return EXIT_FAILURE;
  }


  std::cout << "<<INFO>> DATA LOADED in " << t << std::endl;
  t.reset();
  std::cout << "<<INFO>> DATA Y dimension :" << vs_data->parametrization().dimY() << std::endl;



  // Load a function file representing the BRDF
  ptr<function> brdf = ptr<function>(plugins_manager::load_function(args["brdf"]));
  if(brdf.get() == NULL)
  {  
    std::cout << "<<ERROR>> Could not load the BRDF. Check your file. ABORTING !!!! "  << std::endl;
    return EXIT_FAILURE;  
  }
  std::cout << "<<INFO>> BRDF File Loaded. Starting to compute statistics ... " << std::endl;


  //Conversion 
  ptr<data> generic_data = dynamic_pointer_cast<data>( vs_data );

  vertical_segment*   conv_vs = NULL;
  bool conversion_necessary = true;
  std::cout << "<<INFO>> Converting data to function parametrization if needed" << std::endl;

  
  t.start();
  convertDataToFunctionParam( generic_data, brdf, conversion_necessary, conv_vs );
  t.stop();
  std::cout << "<<INFO>> Data converted in " << t << std::endl;
  t.reset();

  ptr<data> converted_data = NULL;
  if( conversion_necessary )
  {
    ptr<vertical_segment> sp_conv_vs( conv_vs );
    converted_data = dynamic_pointer_cast<data>(sp_conv_vs ); 
  }
  else
  {
    converted_data = generic_data;
  }

  //Here we go new way and faster because we evaluate the function just once
  Eigen::ArrayXXd data_x =  Eigen::ArrayXXd::Zero( converted_data->size(), converted_data->parametrization().dimX() );
  Eigen::ArrayXXd data_y =  Eigen::ArrayXXd::Zero( converted_data->size(), converted_data->parametrization().dimY() ) ;

  t.start();
  demultiplexData( converted_data, data_x, data_y );
  t.stop();
  std::cout << "<<INFO>> Demultiplex data in " << t << std::endl;
  t.reset();

  Eigen::ArrayXXd f_y    = Eigen::ArrayXXd::Zero( converted_data->size(), brdf->dimY() ) ;

  t.start();
  evaluateFunctionAtData( data_x, brdf, f_y );
  t.stop();
  std::cout << "<<INFO>> Evaluate function for all data in  " << t << std::endl;
  t.reset();
   
  Eigen::VectorXd L1_norm   = Eigen::VectorXd::Zero( brdf->dimY() );
  Eigen::VectorXd L2_norm   = Eigen::VectorXd::Zero( brdf->dimY() );
  Eigen::VectorXd L3_norm   = Eigen::VectorXd::Zero( brdf->dimY() );
  Eigen::VectorXd LInf_norm = Eigen::VectorXd::Zero( brdf->dimY() );

  Eigen::VectorXd mse       = Eigen::VectorXd::Zero( brdf->dimY() );
  Eigen::VectorXd rmse      = Eigen::VectorXd::Zero( brdf->dimY() );

  t.start();
  fastNormComputation( data_y, f_y, L1_norm, L2_norm, L3_norm, LInf_norm, mse, rmse);
  t.stop();
  std::cout << "<<INFO>> Fast Norm Computations in  " << t << std::endl;
  std::cout << "<<INFO>> L1_norm " << L1_norm << std::endl
            << "<<INFO>> L2_norm " << L2_norm << std::endl
            << "<<INFO>> L3_norm " << L3_norm << std::endl
            << "<<INFO>> Linf_norm " << LInf_norm << std::endl
            << "<<INFO>> Mse  " << mse << std::endl
            << "<<INFO>> Rmse " << rmse << std::endl;
  t.reset();




  
  //Norm L1
  // t.start(); 
  // vec const  norm_l1 = Norm::L1(converted_data, brdf);
  // t.stop();
  // std::cout << "<<INFO>> From Norm L1:" << norm_l1 << " (computed in " << t << ")" << std::endl;  
  // t.reset();

  //Norm L2
  // t.start();
  // vec const norm_l2 = Norm::L2(converted_data, brdf);
  // t.stop();
  // std::cout << "<<INFO>> From Norm L2:" << norm_l2 << " (computed in " << t << ")" << std::endl;
  // t.reset();

  //Norm L3
  // t.start();
  // vec const norm_l3 = Norm::Lp(converted_data, brdf, 3.0);
  // t.stop();
  // std::cout << "<<INFO>> From Norm L3:"    << norm_l3 << " (computed in " << t << ")" << std::endl;
  // t.reset();

  //Norm LInf
  // t.start();
  // vec const norm_linf = Norm::LInf(converted_data, brdf);
  // t.stop();
  // std::cout << "<<INFO>> From Norm LINF: " <<  norm_linf << " (computed in " << t << ")" << std::endl;
  // t.reset();

  //MSE and RMSE
  // t.start();
  // vec const mse = Error::meanSquareError(converted_data, brdf);
  // t.stop();
  // std::cout << "<<INFO>> Mean Square Error = " << mse << " (computed in " << t << ")" << std::endl;
  // t.reset();

  // t.start();
  // vec const rmse = Error::rootMeanSquareError(converted_data, brdf);
  // t.stop();
  // std::cout << "<<INFO>> Root Mean Square Error = " << rmse << " (computed in " << t << ")" << std::endl;
  // t.reset();

  // Compute Weighted Norms by cosine factors
  vec cosine_theta_view;
  vec cosine_theta_light;
  computeCosineFactorsFromData( converted_data, cosine_theta_light, cosine_theta_view );
  vec cosine_light_view = cosine_theta_light.cwiseProduct( cosine_theta_view );
  


  Eigen::VectorXd w_cosine_light_mse;
  Eigen::VectorXd w_cosine_light_rmse;

  fastWeightedErrors( data_y, f_y, cosine_theta_light, w_cosine_light_mse, w_cosine_light_rmse);

  std::cout << "<<INFO>> Weighted MSE by cosine of the light direction: " << w_cosine_light_mse << std::endl;
  std::cout << "<<INFO>> Weighted Root-MSE by cosine of the light direction: " << w_cosine_light_rmse << std::endl;


  Eigen::VectorXd w_cosine_light_view_mse;
  Eigen::VectorXd w_cosine_light_view_rmse;
  fastWeightedErrors( data_y, f_y, cosine_light_view, w_cosine_light_view_mse, w_cosine_light_view_rmse);

  std::cout << "<<INFO>> Weighted MSE by cos(theta_light) cos(theta_view): " << w_cosine_light_view_mse << std::endl;
  std::cout << "<<INFO>> Weighted Root-MSE by cos(theta_light) cos(theta_view): " << w_cosine_light_view_rmse << std::endl;



  // t.start();
  // vec const norm_l2_cos_light = Norm::weightedL2(converted_data, brdf, cosine_theta_light);
  // t.stop();
  // std::cout << "<<INFO>> Weighted L2 with BRDF*cos(theta_light) " 
  //           << norm_l2_cos_light << " (computed in " << t << ")" << std::endl;
  // t.reset();

  // t.start();
  // vec const norm_l2_cos_light_view = Norm::weightedL2(converted_data, brdf, cosine_light_view);
  // t.stop();
  // std::cout << "<<INFO>> Weighted L2 with BRDF*cos(theta_light)*cos(theta_view) " 
  //           << norm_l2_cos_light_view << " (computed in " << t << ")" << std::endl;
  // t.reset();


  //If output is not void we output the different metrics to a file
  if( args.is_defined("output") )
  {
    std::string output_filename = args["output"];
    if( output_filename != "")
    {
      std::ofstream  fwriter(output_filename.c_str() );

      if( fwriter.good() )
      {
        fwriter << "#FIT2STATS EXPORT. DIMENSIONS ARE THE ONE FROM THE BRDF FUNCTION" << std::endl;
        fwriter << "#CMD " << args.get_cmd() << std::endl;
        fwriter << "#DIM " << brdf->dimX() << " " << brdf->dimY() << std::endl;
        fwriter << std::endl;

        fwriter << "L1 :" << L1_norm << std::endl;
        fwriter << "L2 :" << L2_norm << std::endl;
        fwriter << "L3 :" << L3_norm << std::endl;
        fwriter << "LINF :" << LInf_norm << std::endl;
        fwriter << "MSE :" << mse << std::endl;
        fwriter << "RMSE :" << rmse << std::endl;
        fwriter << "MSE * cos(theta_light) :" << w_cosine_light_mse << std::endl;
        fwriter << "RMSE * cos(theta_light) :" << w_cosine_light_rmse << std::endl;
        fwriter << "MSE * cos(theta_light) * cos(theta_view) :" << w_cosine_light_mse << std::endl;
        fwriter << "RMSE * cos(theta_light) cos(theta_view)  :" << w_cosine_light_rmse << std::endl;
        fwriter << std::endl;
      }
      else
      {
        std::cerr << "<<ERROR>> Could not save metrics to " << output_filename << std::endl;
        return EXIT_FAILURE;
      }
    }
  }


  //Comparisons with the norm methods provided in function.cpp
  t.start();
  double const L2   = brdf->L2_distance( generic_data ) ;
  t.stop();
  std::cout << "<<INFO>> L2, (computed from function)  distance to data = " << L2  
            << " (computed in " << t << ")" << std::endl;
  t.reset();

  t.start();
  double const Linf = brdf->Linf_distance( generic_data );  
  t.stop();
  std::cout << "<<INFO>> Linf distance to data = " << Linf 
            << " (computed in " << t << ")" << std::endl;
  t.reset();

  return EXIT_SUCCESS;
}
