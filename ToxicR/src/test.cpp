#include <immintrin.h>
#include <iostream>
#ifdef R_COMPILATION
//necessary things to run in R    
  #include <RcppEigen.h>
  #include <RcppGSL.h>
#else 
  #include <Eigen/Dense>
#endif

#include <normal_FUNL_NC.h>
#include "cBMDstatmod.h"
#include <IDPrior.h>

#include <gsl/gsl_randist.h>
using namespace Rcpp;
using namespace std; 
// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector test_skn(Eigen::MatrixXd Y, Eigen::MatrixXd D) {
    
    
    Eigen::MatrixXd prior(7,5); 
    prior << 1 ,0 , 5 , -100 , 100,
             1 ,0 , 5 , -100 , 100,
             2 ,0 , 1, 0, 100,
             2 , 0 , 1  , 0   , 10,
             2 , 0, 1,  0, 1.0,
             2 , 0 , 1  , 0   , 10,
             1,  0 , 1    ,-18  , 18;
    
    std::vector<bool> fixedB(prior.rows());
    std::vector<double> fixedV(prior.rows());
    for (int i = 0; i < prior.rows(); i++) {
      fixedB[i] = false;
      fixedV[i] = 0.0;
    }
    // value to return
    IDcontinuousPrior model_prior(prior);
    normalFUNL_BMD_NC likelihood(Y, D,false,true,0);
   
    // create the Continuous BMD modelds
    cBMDModel<normalFUNL_BMD_NC, IDcontinuousPrior>  model(likelihood, model_prior, fixedB, fixedV, true);	
    
    optimizationResult OptRes = findMAP<normalFUNL_BMD_NC, IDcontinuousPrior>(&model);
    cout << "BMRs" << endl; 
    double bmd = likelihood.bmd_hybrid_extra(OptRes.max_parms,0.1,true,0.01); 
    cout << OptRes.max_parms << endl;
    Eigen::MatrixXd result = profile_cBMDNC<normalFUNL_BMD_NC, IDcontinuousPrior>(  &model,
                                                                                    CONTINUOUS_BMD_HYBRID_EXTRA,
                                                                                     bmd,0.10,0.01,
                                                                                     0.02,
                                                                                     1.920729+0.1,  // the plus 0.1 is to go slightly beyond
                                                                                     true);
    cout << endl << bmd << endl; 
    //Rcpp::Rcout << result << std::end; 
      
          
    /*cout << likelihood.bmd_point(OptRes.max_parms,1.8,false) << endl; 
    cout << likelihood.bmd_stdev(OptRes.max_parms,2.2,false) << endl; 
    cout << likelihood.bmd_hybrid_extra(OptRes.max_parms, 0.05, false , 0.01) << endl; 
    */
    return wrap(result); 
}


