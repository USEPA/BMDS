

#include <iostream>
#include <string>
#include <vector>

//#include "bmds_entry.h"
//#include "continuous_model_functions.h"


//const Map<MatrixXd> A(as<Map<MatrixXd>>(AA));

//#include  <statmod.h>


//#include <normal_likelihoods.h>
//#include <normalModels.h>
#define STRICT_R_HEADERS

//#include <statmod.h>
//#include <dBMDstatmod.h>
//#include <log_likelihoods.h>
//#include <binomModels.h>
//#include <IDPrior.h>

//#include <bmd_calculate.h>


#include <DichHillBMD_NC.h>
#include <DichMultistageBMD_NC.h>
#include <DichLogLogisticBMD_NC.h>
#include <DichLogProbitBMD_NC.h>
#include <DichWeibullBMD_NC.h>
#include <DichGammaBMD_NC.h>
#include <DichQlinearBMD_NC.h>
#include <DichLogisticBMD_NC.h>
#include <DichProbitBMD_NC.h>

#include "mcmc_analysis.h"

#include <bmd_calculate.h>

#include "normal_HILL_NC.h"
#include "normal_POWER_NC.h"
#include "normal_POLYNOMIAL_NC.h"
#include "normal_EXP_NC.h"

#include "lognormal_HILL_NC.h"
#include "lognormal_POWER_NC.h"
#include "lognormal_POLYNOMIAL_NC.h"
#include "lognormal_EXP_NC.h"

#include "continuous_clean_aux.h"

#ifdef R_COMPILATION
    //necessary things to run in R
    #include <RcppEigen.h>
    #include <RcppGSL.h>
#else
    #include <Eigen/Dense>
#endif

#include "bmds_entry.h"
#include "bmdStruct.h"
#include "continuous_clean_aux.h"
#include "continuous_entry_code.h"
#include "dichotomous_entry_code.h"
    
#include "list_r_conversion.h"

using namespace Rcpp;
using namespace std;
using Eigen::Map;
using Eigen::MatrixXd;
using Rcpp::as;

Eigen::MatrixXd fix_sample(Eigen::MatrixXd A, dich_model mtype, double max){

  // Note: Samples are by column. 
  switch (mtype){
    case dich_model::d_hill:
      A.row(2).array() +=  A.row(3).array()*log(1./max); 
      break;
    case dich_model::d_gamma:
      A.row(2) *= 1./max; 
      break;
    case dich_model::d_logistic:
      A.row(1) *= 1./max; 
      break;
    case dich_model::d_loglogistic:
      A.row(1).array() += A.row(2).array()*log(1./max);  
      break;
    case dich_model::d_logprobit:
      A.row(1).array() +=  A.row(2).array()*log(1./max); 
      break;
    case dich_model::d_multistage:
        for (int j = 1; j < A.rows(); j++){
          A.row(j) *= pow(1./max,j);
        }
      break;
    case dich_model::d_probit:
      A.row(1) *= 1./max; 
      break;
    case dich_model::d_qlinear:
      A.row(1) *= 1./max; 
      break;
    default:
      for( int i = 0; i < A.cols(); i++){
            A(2,i) *= pow( 1/max, A(1,i)); 
      }
  }
  
  return A;   
}

void bmd_single_dichotomous_mcmc_fitter(const dichotomous_analysis *input,
                                        bmd_analysis_MCMC * output)
{

  
  Eigen::MatrixXd Y(input->n,2); 
  Eigen::MatrixXd D(input->n,1);
  Eigen::MatrixXd pr(input->parms,5);
  std::vector<bool>   fixedB(pr.rows());
  std::vector<double> fixedV(pr.rows()); 
  
  double max_dose = input->doses[0]; 

  for (int i = 0; i < input->n; i++){
    Y(i,0) = input->Y[i];
    Y(i,1) = input->n_group[i]; 
    D(i,0) = input->doses[i]; 
    if (input->doses[i] > max_dose){
      max_dose = input->doses[i]; 
    }
  }

  D = (1/max_dose)*D; 
  for (int i = 0; i < input->parms; i++){
    for (int j = 0; j < 5; j++){
      pr(i,j) = input->prior[i + j*input->parms]; 
    }
  }
  
  
  int degree = 1;
  for(int i = 0; i < pr.rows(); i++){
    fixedB[i] = false;
    fixedV[i] = 0;
  }
  mcmcSamples Rval;

  switch (input->model ){
  case dich_model::d_hill:
    Rval = MCMC_bmd_analysis_DNC< dich_hillModelNC,IDcontinuousPrior> (Y, D, pr, 
                                                                       fixedB,  fixedV, degree,
                                                                       input->BMR, input->BMD_type, 
                                                                       input->alpha,input->samples);
    break;
  case dich_model::d_gamma:
    Rval = MCMC_bmd_analysis_DNC< dich_gammaModelNC,IDcontinuousPrior> (Y, D, pr,
                                                                        fixedB,  fixedV, degree,
                                                                        input->BMR, input->BMD_type, 
                                                                        input->alpha,input->samples);
    break;
  case dich_model::d_logistic:
    Rval = MCMC_bmd_analysis_DNC< dich_logisticModelNC,IDcontinuousPrior>(Y, D, pr,
                                                                          fixedB,  fixedV, degree,
                                                                          input->BMR, input->BMD_type, 
                                                                          input->alpha,input->samples);
    
    break;
  case dich_model::d_loglogistic:
    Rval = MCMC_bmd_analysis_DNC< dich_loglogisticModelNC,IDcontinuousPrior> (Y, D, pr,
                                                                              fixedB,  fixedV, degree,
                                                                              input->BMR, input->BMD_type, 
                                                                              input->alpha,input->samples);
    break;
  case dich_model::d_logprobit:
    Rval = MCMC_bmd_analysis_DNC<dich_logProbitModelNC,IDcontinuousPrior> (Y, D, pr,
                                                                           fixedB,  fixedV, degree,
                                                                           input->BMR, input->BMD_type, 
                                                                           input->alpha,input->samples);
    break;
  case dich_model::d_multistage:
    degree = pr.rows() -1; 
    Rval = MCMC_bmd_analysis_DNC< dich_multistageNC,IDcontinuousPrior> (Y, D, pr,
                                                                        fixedB,  fixedV, degree,
                                                                        input->BMR, input->BMD_type, 
                                                                        input->alpha,input->samples);
    break;
  case dich_model::d_probit:
    Rval = MCMC_bmd_analysis_DNC<dich_probitModelNC,IDcontinuousPrior> (Y, D, pr,
                                                                        fixedB,  fixedV, degree,
                                                                        input->BMR, input->BMD_type, 
                                                                        input->alpha,input->samples);
    break;
  case dich_model::d_qlinear:
    Rval =  MCMC_bmd_analysis_DNC<dich_qlinearModelNC,IDcontinuousPrior> (Y, D, pr,
                                                                          fixedB,  fixedV, degree,
                                                                          input->BMR, input->BMD_type, 
                                                                          input->alpha,input->samples);
    break;
  default:
    Rval = MCMC_bmd_analysis_DNC< dich_weibullModelNC,IDcontinuousPrior> (Y, D, pr,
                                                                          fixedB,  fixedV, degree,
                                                                          input->BMR, input->BMD_type, 
                                                                          input->alpha,input->samples);
  break;
  }
  
  Rval.samples =  fix_sample(Rval.samples, input->model, max_dose); 
  output->samples = input->samples; 
  output->model = input->model; 
  for (int i = 0; i < input->samples; i++){
    output->BMDS[i] = Rval.BMD(0,i)*max_dose; // rescale 
    for (int j = 0; j < input->parms;j++){
      output->parms[j+i*input->parms] = Rval.samples(j,i); 
    }
  }
 
  return; 
    
}




//void bmd_single_continuous_mcmc_fitter()

// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::depends(RcppEigen)]]
//////////////////////////////////////////////////////////////////////////
// function: run_dichotomous_single_mcmc
// purpose: takes input, which is assumed to be correct (i.e., filtered
// correctly by the R calling function), and then calls the library to
// run the corresponding analysis. Does MCMC sample
// output: BMD analysis with the model specified by NumericVector model
// [[Rcpp::export]]
List run_dichotomous_single_mcmc(NumericVector model,
              				     Eigen::MatrixXd Y, Eigen::MatrixXd D,
					                 Eigen::MatrixXd pr, NumericVector options){

	dichotomous_analysis mcmcAnal; 
	mcmcAnal.BMD_type =  eExtraRisk;// (options[0]==1)?eExtraRisk:eAddedRisk;
	mcmcAnal.BMR      = options[0]; 
	mcmcAnal.alpha     = options[1];
	mcmcAnal.samples  = options[2]; 
	mcmcAnal.burnin   = options[3]; 
  mcmcAnal.parms    = pr.rows(); 
  mcmcAnal.model = (dich_model)model[0]; 
  mcmcAnal.Y       = new double[Y.rows()] ; 
  mcmcAnal.n_group = new double[Y.rows()] ; 
  mcmcAnal.doses   = new double[D.rows()] ; 
  mcmcAnal.prior   = new double[pr.cols()*pr.rows()];
  mcmcAnal.prior_cols = pr.cols(); 
  mcmcAnal.n          = Y.rows(); 
  mcmcAnal.degree = 0; 
  
  if (mcmcAnal.model == dich_model::d_multistage){
    mcmcAnal.degree = mcmcAnal.parms - 1; 
    cerr << mcmcAnal.degree << endl; 
  }
  
  bmd_analysis_MCMC output; 
  output.samples = mcmcAnal.samples; // initialize
  output.model = (dich_model)0; 
  output.BMDS =  new double[mcmcAnal.samples]; 
  output.parms = new double[mcmcAnal.samples*pr.rows()]; 
  
  for (int i = 0; i < Y.rows(); i++){
    mcmcAnal.Y[i] = Y(i,0); 
    mcmcAnal.n_group[i] = Y(i,1); 
  }
  
  for (int i = 0; i < D.rows(); i++){
    mcmcAnal.doses[i] = D(i,0); 
  }

  // copy in column major order  
  for (int i = 0; i < pr.rows(); i++){
    for (int j = 0; j < pr.cols(); j++){
      mcmcAnal.prior[i + j*pr.rows()] = pr(i,j); 
    }
  }


  dichotomous_model_result res; 
  res.parms = new double[pr.rows()]; 
  res.cov   = new double[pr.rows()*pr.rows()]; 
  res.dist_numE = 200; 
  res.bmd_dist = new double[res.dist_numE*2]; 
 
  estimate_sm_mcmc(&mcmcAnal, &res, &output); 
  
  List rV = convert_dichotomous_fit_to_list(&res); 
  List t2 = convert_MCMC_fit_to_list(&output);
  
  List data_out  = List::create(Named("mcmc_result")=t2,
                                Named("fitted_model")=rV); 
  
  delete(output.BMDS); 
  delete(output.parms); 
  delete(mcmcAnal.Y); 
  delete(mcmcAnal.n_group); 
  delete(mcmcAnal.doses); 
  delete(mcmcAnal.prior); 
  delete(res.parms);   
  delete(res.cov);    
  delete(res.bmd_dist);  
  return data_out;
}



//////////////////////////////////////////////////////////////////////////
// function: run_dichotomous_single_mcmc
// purpose: takes input, which is assumed to be correct (i.e., filtered
// correctly by the R calling function), and then calls the library to
// run the corresponding analysis. Does MCMC sample
// output: BMD analysis with the model specified by NumericVector model

// [[Rcpp::export]]
List run_continuous_single_mcmc(NumericVector model,
                                Eigen::MatrixXd Y, Eigen::MatrixXd D,
                                Eigen::MatrixXd priors, NumericVector options,
                                bool is_logNormal,bool suff_stat){
  unsigned int samples = (unsigned int) options[7]; 
  unsigned int burnin  = (unsigned int) options[8];
  double tail_p = (double) options[6]; 
  bool bConstVar = (bool)options[5]; // check if it is constant variance
  bool is_increasing = (bool)options[4];
  double alpha = (double)options[3];
  double bk_prob = (double)options[2];
  double bmrf  = (double)options[1];
  int riskType = (int)options[0];
  
  continuous_analysis *mcmcAnal = new continuous_analysis; 
  distribution dtype; 
  if (is_logNormal){
    dtype = distribution::log_normal; 
  }else{
    if (bConstVar){
      dtype = distribution::normal;
    }else{
      dtype = distribution::normal_ncv; 
    }
  }
  
  mcmcAnal->model   = (cont_model) model[0]; 
  mcmcAnal->Y       =    new double[Y.rows()]; 
  mcmcAnal->n       =    Y.rows(); 
  mcmcAnal->n_group =    new double[Y.rows()]; 
  mcmcAnal->sd      =    new double[Y.rows()]; 
  mcmcAnal->doses   =    new double[Y.rows()]; 
  mcmcAnal->prior   =    new double[priors.rows()*priors.cols()]; 
  mcmcAnal->isIncreasing = is_increasing; 
  mcmcAnal->disttype     = dtype; 
  mcmcAnal->prior_cols   = priors.cols(); 
  mcmcAnal->parms        = priors.rows(); 
  mcmcAnal->alpha        = alpha; 
  mcmcAnal->BMD_type     = riskType; 
  mcmcAnal->BMR          = bmrf; 
  mcmcAnal->samples      = samples; 
  mcmcAnal->burnin       = burnin; 
  mcmcAnal->tail_prob    = tail_p; 
  mcmcAnal->suff_stat    = suff_stat; 
  
  bmd_analysis_MCMC  *output = new bmd_analysis_MCMC; 
  output->parms = new double[samples*mcmcAnal->parms]; 
  output->BMDS  = new double[samples]; 
  
  ///////////////////////////////////////////////////////////////////////
  if (suff_stat){
    Eigen::MatrixXd tY = Y.col(1);
    Y.col(1) = Y.col(2);
    Y.col(2) = tY;
  }
  
  for (int i = 0; i < Y.rows(); i++){
    mcmcAnal->Y[i] = Y(i,0); 
    mcmcAnal->doses[i] = D(i,0); 
    if (suff_stat){
      mcmcAnal->n_group[i] = Y(i,2);
      mcmcAnal->sd[i]      = Y(i,1); 
    }
  }

  for (int i = 0; i < mcmcAnal->parms; i++ ){
    for (int j = 0; j < mcmcAnal->prior_cols; j++){
      mcmcAnal->prior[i+j*mcmcAnal->parms] = priors(i,j); 
    }
  }
  ////////////////////////////////////
  continuous_model_result *res = new_continuous_model_result( mcmcAnal->model,
                                                              mcmcAnal->parms,
                                                              200);
  
  estimate_sm_mcmc(mcmcAnal,
                   res     ,
                   output) ;
  
  
  List rV = convert_continuous_fit_to_list(res); 
  List t2 = convert_MCMC_fit_to_list(output);
  List data_out  = List::create(Named("mcmc_result")=t2,
                                Named("fitted_model")=rV); 

  del_mcmc_analysis(output);
  del_continuous_model_result(res); 
  del_continuous_analysis(*mcmcAnal);

  return wrap(data_out);
}


 
