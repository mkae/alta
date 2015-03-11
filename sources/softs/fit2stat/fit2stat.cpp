/* ALTA --- Analysis of Bidirectional Reflectance Distribution Functions

   Copyright (C) 2013, 2014 Inria

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

// TODO: Finish the implementation
// TODO : Write the documentation
// TODO: Move this class to ALTA Core
class Norm
{
  public:

  static vec L2( ptr<data> const & data, ptr<function> const & f)
  {
    vec l2 = vec::Zero( f->dimY() );

    if( data->input_parametrization() == f->input_parametrization() )
    {
      if( data->output_parametrization() == f->output_parametrization() )
      {
        //All parametrization are the same. Compute the Metric        
        for( unsigned int i=0; i < data->size(); i++)
        {
          vec const dat = data->get(i);          
          vec const data_x = dat.head( data->dimX() );
          vec const data_y = dat.tail( data->dimY() );

          l2 += (f->value(data_x) - data_y).array().square().matrix();
        }

        l2 /= static_cast<double>( data->size() );
        return l2;
      }
      else //Output parametrizations are not the same
      {
        assert(0);
        NOT_IMPLEMENTED();
        return l2;
      }
    }
    else // Input Param are not the same
    {
      if( data->output_parametrization() == f->output_parametrization() )
      {
        // Only input parametrizations are different
        // Convert data input parametrization to function input parametrization
        for( unsigned int i=0; i < data->size(); i++)
        {
          vec const dat = data->get(i);          
          vec const data_y = dat.tail( data->dimY() );

          vec c_data_x = vec::Zero( f->dimX() );

          params::convert(&dat[0], data->input_parametrization(), f->input_parametrization(), &c_data_x[0] );
  
          l2 += (f->value(c_data_x) - data_y).array().square().matrix();
        }

        l2 /= static_cast<double>( data->size() );
        return l2;
      }
      else
      {
        //All parametrizations are different
        assert(0);
        NOT_IMPLEMENTED();
        return l2;
      }      
    }//End of else on input parametrization
  }

  static double weightedL2( ptr<data> const & data,  ptr<function> const & f, vec const & weights )
  {
    //Check the size of the weights are equal
    assert( data->size() == weights.size()  );
    NOT_IMPLEMENTED();

    //double const weights_sum = weights.sum();

    double weighted_l2 = 0.0;
    return weighted_l2;
  }

  static double LInf()
  {
    assert(0);
    NOT_IMPLEMENTED();

    double linf = 0.0;
    return linf;
  }

  static vec squareRootL2( ptr<data> const & data, ptr<function> const & f )
  {
    vec res = L2( data, f );
    return res.cwiseSqrt();
  }
  
  static vec cubicRootL2( ptr<data> const & data, ptr<function> const & f)
  {
    vec res = L2( data, f );
    return res.array().pow(1.0 / 3.0);
  }

};




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

    return EXIT_FAILURE;
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
  ptr<data> vs_data = new vertical_segment();
  //ptr<data> vs_data = plugins_manager::get_data("vertical_segment");
  vs_data->load(args["input"], args);
  
  std::cout << "<<INFO>> DATA LOADED " << std::endl;



  // Load a function file
  ptr<function> brdf = plugins_manager::get_function(args["brdf"]) ;
  if(brdf.get() == NULL)  {  return EXIT_FAILURE;  }
  std::cout << " BRDF File Loaded. Starting to compute statistics ... " << std::endl;

  // Check the compatibility between the data and the function
  // Quite useless now. Will be removed
  plugins_manager::check_compatibility(vs_data, brdf, args);
  
  
  double const L2   = brdf->L2_distance(vs_data);
  std::cout << "<<INFO>> L2   distance to data = " << L2   << std::endl;
  
  std::cout << "<<INFO>> From NORM: L2   = " << Norm::L2(vs_data, brdf)   << std::endl;
  std::cout << "<<INFO>> From NORM: sqrt(L2) = " << Norm::squareRootL2(vs_data, brdf)   << std::endl;
  std::cout << "<<INFO>> From NORM: L2^(1/3) = " << Norm::cubicRootL2(vs_data, brdf) << std::endl;


  double const Linf = brdf->Linf_distance(vs_data);  
  std::cout << "<<INFO>> Linf distance to data = " << Linf << std::endl;





  return EXIT_SUCCESS;
}