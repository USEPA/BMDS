#ifdef R_COMPILATION  
#include <RcppEigen.h>
  #include <RcppGSL.h>
#else 
  #include <Eigen/Dense>
#endif

#include <gsl/gsl_randist.h>
#include "bmdStruct.h"
#include <math.h>
#include <stdio.h>

#include <string>
#include <vector>
#include <algorithm>
#include <limits>

#include "DichHillBMD_NC.h"
#include "DichMultistageBMD_NC.h"
#include "DichLogLogisticBMD_NC.h"
#include "DichLogProbitBMD_NC.h"
#include "DichWeibullBMD_NC.h"
#include "DichGammaBMD_NC.h"
#include "DichQlinearBMD_NC.h"
#include "DichLogisticBMD_NC.h"
#include "DichProbitBMD_NC.h"
#include "IDPrior.h"

#include "bmds_entry.h"
#include "bmdStruct.h"

#include "bmd_calculate.h"
#include "dichotomous_entry_code.h"
#include "mcmc_analysis.h"



void transfer_dichotomous_model(bmd_analysis a, dichotomous_model_result *model){
  if (model){
    model->nparms = a.COV.rows(); 
    model->max = a.MAP; 
    for (int i = 0; i< model->dist_numE; i ++){
      double temp = double(i)/double(model->dist_numE); 
      model->bmd_dist[i] = a.BMD_CDF.inv(temp);     // BMD @ probability
     // cerr << model->bmd_dist[i] << ":" << temp << endl; 
      model->bmd_dist[model->dist_numE + i] = temp; // probability 
    }
    
    for (int i = 0; i < model->nparms; i++){
      model->parms[i] = a.MAP_ESTIMATE(i,0); 
      for (int j = 0; j < model->nparms; j++){
        model->cov[i + j*model->nparms] = a.COV(i,j); 
      }
    }
  }
}


void estimate_ma_MCMC(dichotomousMA_analysis *MA,
                      dichotomous_analysis   *DA,
                      dichotomousMA_result   *res,
                      ma_MCMCfits            *ma){ 
  return ; 
}

void estimate_ma_laplace(dichotomousMA_analysis *MA,
                         dichotomous_analysis *DA ,
                         dichotomousMA_result *res){
  return ; 
}

void estimate_sm_laplace(dichotomous_analysis *DA ,
                         dichotomous_model_result *res){
  return ; 
}

void estimate_sm_mcmc(dichotomous_analysis *DA, 
                      dichotomous_model_result *res,
                      bmd_analysis_MCMC *mcmc){
 
  ///////////////////////////////////
  Eigen::MatrixXd Y(DA->n,2); 
  Eigen::MatrixXd D(DA->n,1); 
  Eigen::MatrixXd prior(DA->parms,DA->prior_cols);
  for (int i = 0; i < DA->n; i++){
      Y(i,0) = DA->Y[i]; Y(i,1) = DA->n_group[i]; 
      D(i,0) = DA->doses[i]; 
  }
  for (int i = 0; i < DA->parms; i++){
      for (int j = 0; j < DA->prior_cols; j++){
        prior(i,j) = DA->prior[i + j*DA->parms]; 
      } 
  }  // copy the prior over. 
  mcmcSamples a; 
  std::vector<bool> fixedB; 
  std::vector<double> fixedV; 
  for (int i = 0; i < prior.rows(); i++){
    fixedB.push_back(false);
    fixedV.push_back(0.0); 
  }
  switch (DA->model){
    case dich_model::d_hill:
      a =  MCMC_bmd_analysis_DNC<dich_hillModelNC,IDPrior> (Y,D,prior,
                                     fixedB, fixedV, DA->degree,
                                     DA->BMR, DA->BMD_type, DA->alpha, DA->samples);
    break; 
    case dich_model::d_gamma:
   
      a =  MCMC_bmd_analysis_DNC<dich_gammaModelNC,IDPrior> (Y,D,prior,
                                                            fixedB, fixedV, DA->degree,
                                                            DA->BMR, DA->BMD_type, DA->alpha, DA->samples);
  
    break; 
    case dich_model::d_logistic:
      a =  MCMC_bmd_analysis_DNC<dich_logisticModelNC,IDPrior> (Y,D,prior,
                                                            fixedB, fixedV, DA->degree,
                                                            DA->BMR, DA->BMD_type, DA->alpha, DA->samples);
    break; 
    case dich_model::d_loglogistic:
      a =  MCMC_bmd_analysis_DNC<dich_loglogisticModelNC,IDPrior> (Y,D,prior,
                                                            fixedB, fixedV, DA->degree,
                                                            DA->BMR, DA->BMD_type, DA->alpha, DA->samples);
    break;
    case dich_model::d_logprobit:
      a =  MCMC_bmd_analysis_DNC<dich_logProbitModelNC,IDPrior> (Y,D,prior,
                                                            fixedB, fixedV, DA->degree,
                                                            DA->BMR, DA->BMD_type, DA->alpha, DA->samples);
    break; 
    case dich_model::d_multistage:
      a =  MCMC_bmd_analysis_DNC<dich_multistageNC,IDPrior> (Y,D,prior,
                                                            fixedB, fixedV, DA->degree,
                                                            DA->BMR, DA->BMD_type, DA->alpha, DA->samples);
    break;
    case dich_model::d_probit: 
      a =  MCMC_bmd_analysis_DNC<dich_probitModelNC,IDPrior> (Y,D,prior,
                                                            fixedB, fixedV, DA->degree,
                                                            DA->BMR, DA->BMD_type, DA->alpha, DA->samples);
    break;
    case dich_model::d_qlinear: 
      a =  MCMC_bmd_analysis_DNC<dich_qlinearModelNC,IDPrior> (Y,D,prior,
                                                            fixedB, fixedV, DA->degree,
                                                            DA->BMR, DA->BMD_type, DA->alpha, DA->samples);
    break; 
    case dich_model::d_weibull:
      a =  MCMC_bmd_analysis_DNC<dich_weibullModelNC,IDPrior> (Y,D,prior,
                                                             fixedB, fixedV, DA->degree,
                                                             DA->BMR, DA->BMD_type, DA->alpha, DA->samples);
    default: 
    break; 
  }
  bmd_analysis b; 
  b = create_bmd_analysis_from_mcmc(DA->burnin,a);
  mcmc->model = DA->model; 
  mcmc->burnin = DA->burnin; 
  mcmc->samples = DA->samples; 
  mcmc->nparms = DA->parms; 
  transfer_mcmc_output(a,mcmc); 
  res->model = DA->model; 
  transfer_dichotomous_model(b,res);
 
  return; 
}