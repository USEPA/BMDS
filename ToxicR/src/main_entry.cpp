#include <RcppEigen.h>
#include <RcppGSL.h>

#include <iostream>
#include <string>
#include <vector>

#include "bmds_entry.h"


#include  <statmod.h>

#include <log_likelihoods.h>
#include <normal_likelihoods.h>
#include <normalModels.h>
#include <binomModels.h>
#include <IDPrior.h>

#include "bmd_calculate.h"

#include "normal_HILL_NC.h"
#include "normal_POWER_NC.h"
#include "normal_POLYNOMIAL_NC.h"
#include "normal_EXP_NC.h"

#include "lognormal_HILL_NC.h"
#include "lognormal_POWER_NC.h"
#include "lognormal_POLYNOMIAL_NC.h"
#include "lognormal_EXP_NC.h"

#include "continuous_clean_aux.h"
#include "continuous_entry_code.h"
#include "dichotomous_entry_code.h"
#include "list_r_conversion.h"

using namespace Rcpp;
using namespace std;

#define MAX_PARMS 32 // Should never get close to this many!!!

// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::depends(RcppEigen)]]
//////////////////////////////////////////////////////////////////////////
// function: run_single_dichotomous
// purpose: takes input, which is assumed to be correct (i.e., filtered
// correctly by the R calling function), and then calls the library to
// run the corresponding analysis.
// output: BMD analysis with the model specified by NumericVector model
//
// [[Rcpp::export]]
List run_single_dichotomous(NumericVector model,
                            Eigen::MatrixXd data, Eigen::MatrixXd pr,
                            NumericVector options1, IntegerVector options2) 
{

  dichotomous_analysis Anal; 
  Anal.BMD_type =  (options1[0]==1)?eExtraRisk:eAddedRisk;
  Anal.BMR      =  options1[0]; 
  Anal.alpha    = options1[1];
  Anal.parms    = pr.rows(); 
  Anal.model    = (dich_model)model[0]; 
  Anal.Y        = new double[data.rows()] ; 
  Anal.n_group  = new double[data.rows()] ; 
  Anal.doses    = new double[data.rows()] ; 
  Anal.prior    = new double[pr.cols()*pr.rows()];
  Anal.prior_cols = pr.cols(); 
  Anal.n          = data.rows(); 
  Anal.degree = options2[1]; 

  if (Anal.model == dich_model::d_multistage){
    Anal.degree = Anal.parms - 1; 
  }

  for (int i = 0; i < data.rows(); i++){
    Anal.Y[i] = data(i,1); 
    Anal.n_group[i] = data(i,2); 
  }
  
  for (int i = 0; i < data.rows(); i++){
    Anal.doses[i] = data(i,0); 
  }
  
  // copy in column major order  
  for (int i = 0; i < pr.rows(); i++){
    for (int j = 0; j < pr.cols(); j++){
      Anal.prior[i + j*pr.rows()] = pr(i,j); 
    }
  }
  
  dichotomous_model_result res; 
  res.parms = new double[pr.rows()]; 
  res.cov   = new double[pr.rows()*pr.rows()]; 
  res.dist_numE = 200; 
  res.bmd_dist = new double[res.dist_numE*2]; 
  
  estimate_sm_laplace(&Anal, &res); 
  List rV = convert_dichotomous_fit_to_list(&res); 
  
  delete(Anal.Y); 
  delete(Anal.n_group); 
  delete(Anal.doses); 
  delete(Anal.prior); 
  delete(res.parms);   
  delete(res.cov);    
  delete(res.bmd_dist);  
  return rV;
}

// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::depends(RcppEigen)]]
//////////////////////////////////////////////////////////////////////////
// function: run_single_continuous
// purpose: takes input, which is assumed to be correct (i.e., filtered
// correctly by the R calling function), and then calls the library to
// run the corresponding analysis.
// output: BMD analysis with the model specified by NumericVector model
// [[Rcpp::export]]
List run_continuous_single(IntegerVector model, 
                           Eigen::MatrixXd Y, Eigen::MatrixXd X,
                           Eigen::MatrixXd prior, NumericVector options,
                           IntegerVector dist_type){
    
    bool   is_increasing = (bool)options[4]; 	double alpha = (double)options[3];
    double tail_p = (double)options[2]; 	double bmrf  = (double)options[1];
    int    riskType = (int)options[0];   
    unsigned int samples = (unsigned int) options[5];
    
    ////////////////////////////////////////////////
    /// Set up the analysis
    ////////////////////////////////////////////////
    continuous_analysis anal; 
    anal.Y       =    new double[Y.rows()]; 
    anal.n       =    Y.rows(); 
    anal.n_group =    new double[Y.rows()]; 
    anal.sd      =    new double[Y.rows()]; 
    anal.doses   =    new double[Y.rows()]; 
    anal.model   =    (cont_model) model[0]; 
    anal.disttype     = dist_type[0]; 
    anal.isIncreasing = is_increasing; 
    anal.alpha        = 0.005; //alpha for analyses; 
    anal.BMD_type     = riskType; 
    anal.BMR          = bmrf; 
    anal.samples      = samples; 
    anal.tail_prob    = tail_p; 
    anal.suff_stat    = Y.cols()==3;
    anal.parms        = prior.rows();
    anal.prior_cols   = prior.cols(); 
    anal.prior   = new double[prior.rows()*prior.cols()]; 
    cp_prior(prior,anal.prior);
    
    for (int i = 0; i < Y.rows(); i++){
      anal.Y[i] = Y(i,0); 
      anal.doses[i] = X(i,0); 
      if (Y.cols() == 3){ //sufficient statistics
        anal.n_group[i] = Y(i,2);
        anal.sd[i]      = Y(i,1); 
      }
    }
    ////////////////////////////////////
    continuous_model_result *result = new_continuous_model_result( anal.model,
                                                                   anal.parms,
                                                                   200); //have 200 equally spaced values
    ////////////////////////////////////
    estimate_sm_laplace(&anal,result);
    List rV = convert_continuous_fit_to_list(result); 	
    // free up memory
    del_continuous_model_result(result); 
    del_continuous_analysis(anal);
    return rV; 
    
  }
