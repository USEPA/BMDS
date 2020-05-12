

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
	mcmcAnal.BMR     = options[0]; 
	mcmcAnal.alpha   = options[1];
	mcmcAnal.samples = options[2]; 
  mcmcAnal.parms   = pr.rows(); 
  
  
  mcmcAnal.model = (dich_model)model[0]; 
  mcmcAnal.Y       = new double[Y.rows()] ; 
  mcmcAnal.n_group = new double[Y.rows()] ; 
  mcmcAnal.doses   = new double[D.rows()] ; 
  mcmcAnal.prior   = new double[pr.cols()*pr.rows()];
  mcmcAnal.n       = Y.rows(); 
  
  bmd_analysis_MCMC output; 
  output.samples = 0; // initialize
  output.model = (dich_model)0; 
  output.BMDS = new double[mcmcAnal.samples]; 
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
  
  
  
  bmd_single_dichotomous_mcmc_fitter(&mcmcAnal,&output); 
  
  Eigen::MatrixXd BMDS(output.samples,1); 
  Eigen::MatrixXd PARMS(output.samples,pr.rows()); 
  
 // disregard the burnin i.e. i = options[3]
  for (int i = 0; i < output.samples; i++){
    BMDS(i,0) = output.BMDS[i]; 
    for (int j=0; j < pr.rows(); j++){
     PARMS(i,j) = output.parms[j + i* pr.rows()]; 
    }
  }
  
  List data_out  = List::create(Named("BMDS")=BMDS,
                                Named("PARMS")=PARMS); 
  
  
  delete(output.BMDS); 
  delete(output.parms); 
  delete(mcmcAnal.Y); 
  delete(mcmcAnal.n_group); 
  delete(mcmcAnal.doses); 
  delete(mcmcAnal.prior); 
  return data_out;
}



void bmd_single_continuous_mcmc_fitter(const continuous_analysis *input,
                                             bmd_analysis_MCMC * output)
  {
   
 
    int input_cols;
    if (input->suff_stat){
      input_cols = 3; 
    }else{
      input_cols = 1; 
    }
    Eigen::MatrixXd Y(input->n,input_cols); 
    Eigen::MatrixXd D(input->n,1);
    
    
    Eigen::MatrixXd pr(input->parms,input->priorD);
    
    std::vector<bool>   fixedB(pr.rows());
    std::vector<double> fixedV(pr.rows()); 

    double max_dose = input->doses[0]; 
    
    for (int i = 0; i < input->n; i++){
      Y(i,0) = input->Y[i];
      D(i,0) = input->doses[i]; 
      if (input->suff_stat){
        Y(i,1) = input->sd[i]; 
        Y(i,2) = input->n_group[i]; 
      }
      if (input->doses[i] > max_dose){
        max_dose = input->doses[i]; 
      }
    }
   
   
    D = (1/max_dose)*D; 
    for (int i = 0; i < input->parms; i++){
       cerr.flush();
      for (int j = 0; j < input->priorD; j++){
        pr(i,j) = input->prior[i + j*input->parms]; 
      }
 
    }
   
    int degree = 1;
    for(int i = 0; i < pr.rows(); i++){
      fixedB[i] = false;
      fixedV[i] = 0;
    }
    
    mcmcSamples Rval;
    int adverseR = 0; 
    
    switch (input->model ){
    case cont_model::hill:
        if (input->isLognorm){
        
            Rval = MCMC_bmd_analysis_CONTINUOUS_LOGNORMAL< lognormalHILL_BMD_NC,IDcontinuousPrior> (Y, D, pr, 
                                                                                                    fixedB,  fixedV,  input->isIncreasing,                                                                                           input->tail_prob,input->suff_stat,
                                                                                                    input->BMR, input->BMD_type,
                                                                                                    input->alpha,input->samples,0);
        }else{
       
            Rval = MCMC_bmd_analysis_CONTINUOUS_NORMAL< normalHILL_BMD_NC,IDcontinuousPrior> (Y, D, pr, 
                                                                                              fixedB,  fixedV,  input->isIncreasing,
                                                                                              input->tail_prob, input->suff_stat,
                                                                                              input->BMR, input->BMD_type, 
                                                                                              input->const_var,
                                                                                              input->alpha,input->samples,0);
        }
        break;
    case cont_model::exp_3:
      adverseR = input->isIncreasing?NORMAL_EXP3_UP: NORMAL_EXP3_DOWN; 
 
      if (input->isLognorm){
        
        Rval = MCMC_bmd_analysis_CONTINUOUS_LOGNORMAL< lognormalEXPONENTIAL_BMD_NC,IDcontinuousPrior> (Y, D, pr, 
                                                                                                       fixedB,  fixedV,  input->isIncreasing,
                                                                                                       input->tail_prob,input->suff_stat,
                                                                                                       input->BMR, input->BMD_type,
                                                                                                       input->alpha,input->samples,adverseR);
      }else{
   
        Rval = MCMC_bmd_analysis_CONTINUOUS_NORMAL< normalEXPONENTIAL_BMD_NC,IDcontinuousPrior>  (Y, D, pr, 
                                                                                                  fixedB,  fixedV,  input->isIncreasing,
                                                                                                  input->tail_prob, input->suff_stat,
                                                                                                  input->BMR, input->BMD_type, 
                                                                                                  input->const_var,
                                                                                                  input->alpha,input->samples,adverseR);
      }

      break;
    case cont_model::exp_5:
      
      adverseR = input->isIncreasing?NORMAL_EXP5_UP: NORMAL_EXP5_DOWN; 
      if (input->isLognorm){
        
        Rval = MCMC_bmd_analysis_CONTINUOUS_LOGNORMAL< lognormalEXPONENTIAL_BMD_NC,IDcontinuousPrior> (Y, D, pr, 
                                                                                                       fixedB,  fixedV,  input->isIncreasing,
                                                                                                       input->tail_prob,input->suff_stat,
                                                                                                       input->BMR, input->BMD_type,
                                                                                                       input->alpha,input->samples,adverseR);
      }else{
        Rval = MCMC_bmd_analysis_CONTINUOUS_NORMAL< normalEXPONENTIAL_BMD_NC,IDcontinuousPrior>  (Y, D, pr, 
                                                                                                  fixedB,  fixedV,  input->isIncreasing,
                                                                                                  input->tail_prob, input->suff_stat,
                                                                                                  input->BMR, input->BMD_type, 
                                                                                                  input->const_var,
                                                                                                  input->alpha,input->samples,adverseR);
      }
   
      break;
    case cont_model::power:
      if (input->isLognorm){
        
        Rval = MCMC_bmd_analysis_CONTINUOUS_LOGNORMAL< lognormalPOWER_BMD_NC,IDcontinuousPrior> (Y, D, pr, 
                                                                                                 fixedB,  fixedV, input->isIncreasing,
                                                                                                  input->BMR, input->BMD_type, input->tail_prob,
                                                                                                  input->suff_stat,
                                                                                                  input->alpha,input->samples,0);
      }else{
        Rval = MCMC_bmd_analysis_CONTINUOUS_NORMAL< normalPOWER_BMD_NC,IDcontinuousPrior> (Y, D, pr, 
                                                                                           fixedB,  fixedV,  input->isIncreasing,
                                                                                           input->tail_prob, input->suff_stat,
                                                                                            input->BMR, input->BMD_type, 
                                                                                            input->const_var,
                                                                                            input->alpha,input->samples,0);
      }
      break;
      
    
    }
    
    output->samples = input->samples; 
    for (int i = 0; i < input->samples; i++){
       output->BMDS[i] = Rval.BMD(0,i)*max_dose; // rescale 
      for (int j = 0; j < Rval.samples.rows();j++){
        output->parms[j+i*Rval.samples.rows()] = Rval.samples(j,i); 
      }
    }
    
    return; 
    
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
  int samples = (int) options[7]; 
  double tail_p = (double) options[6]; 
  bool bConstVar = (bool)options[5]; // check if it is constant variance
  bool is_increasing = (bool)options[4];
  double alpha = (double)options[3];
  double bk_prob = (double)options[2];
  double bmrf  = (double)options[1];
  int riskType = (int)options[0];
  
  
  bool convert = true; 
  double  max_dose = D.maxCoeff(); 
  
  if (!suff_stat){
    // check to see if it can be converted into sufficient statistics4
    int temp = 0; 
    for (int i = 0; i < D.rows(); i++){
      for (int j = 0 ; j < D.rows(); j++){
        if (D(i,0) == D(j,0)){
          temp++; 
        }
      }
      if( temp == 1){
        convert = false; 
      }
    }
    
    if (convert){
      // we can convert the data
      Eigen::MatrixXd SSTAT = createSuffStat( Y, D, is_logNormal); 
      std::vector<double> uniqueX = unique_list(D);
      Eigen::MatrixXd temp_D(uniqueX.size(),1);
      for (int i = 0; i < uniqueX.size(); i++){
        temp_D(i,0) = uniqueX[i]; 
      }
      D = temp_D; 
      Y = SSTAT; 
      suff_stat = true; 
    }
  }
  
  double divisor = get_diviosor( Y,  D); 
  
  if (suff_stat){
    Y = cleanSuffStat(Y,D,is_logNormal);  
  }else{
    Y = (1/divisor)*Y; // scale the data with the divisor term. 
  }
  
  D = (1/max_dose)*D;
  
  continuous_analysis mcmcAnal; 

  
  mcmcAnal.model   = (cont_model) model[0]; 
  mcmcAnal.Y       =    new double[Y.rows()]; 
  mcmcAnal.n       =    Y.rows(); 
  mcmcAnal.n_group =    new double[Y.rows()]; 
  mcmcAnal.sd      =    new double[Y.rows()]; 
  mcmcAnal.doses   =    new double[Y.rows()]; 
  mcmcAnal.prior   =    new double[priors.rows()*priors.cols()]; 
  mcmcAnal.isIncreasing = is_increasing; 
  mcmcAnal.const_var    = bConstVar; 
  mcmcAnal.isLognorm    = is_logNormal; 
  mcmcAnal.priorD       =  priors.cols(); 
  mcmcAnal.parms        = priors.rows(); 
  mcmcAnal.alpha        = alpha; 
  mcmcAnal.BMD_type     = riskType; 
  mcmcAnal.BMR          = bmrf; 
  mcmcAnal.samples      = samples; 
  mcmcAnal.tail_prob    = tail_p; 
  mcmcAnal.suff_stat    = suff_stat; 
  
  bmd_analysis_MCMC  output; 
  output.parms = new double[samples*mcmcAnal.parms]; 
  output.BMDS = new double[samples]; 
  
  ///////////////////////////////////////////////////////////////////////
  if (suff_stat){
    Eigen::MatrixXd tY = Y.col(1);
    Y.col(1) = Y.col(2);
    Y.col(2) = tY;
  }
  
  for (int i = 0; i < Y.rows(); i++){
    mcmcAnal.Y[i] = Y(i,0); 
    mcmcAnal.doses[i] = D(i,0); 
    if (suff_stat){
      mcmcAnal.n_group[i] = Y(i,2);
      mcmcAnal.sd[i]      = Y(i,1); 
    }
  }

  for (int i = 0; i < mcmcAnal.parms; i++ ){
    for (int j = 0; j < mcmcAnal.priorD; j++){
      mcmcAnal.prior[i+j*mcmcAnal.parms] = priors(i,j); 
    }
  }
  
  
  bmd_single_continuous_mcmc_fitter(&mcmcAnal,
                                    &output);
  
  Eigen::MatrixXd BMDS(output.samples,1); 
  Eigen::MatrixXd PARMS(output.samples,priors.rows()); 
  
  // disregard the burnin i.e. i = options[3]
  for (int i = 0; i < output.samples; i++){
    BMDS(i,0) = output.BMDS[i]*max_dose; 
    for (int j=0; j < priors.rows(); j++){
      PARMS(i,j) = output.parms[j + i* priors.rows()]; 
    }
    Eigen::MatrixXd temp = rescale_parms(PARMS.row(i).transpose(), mcmcAnal.model ,
                                          max_dose, divisor, is_logNormal); 
    PARMS.row(i) = temp.transpose(); 
  }
 
  List data_out  = List::create(Named("BMD")=BMDS,
                                Named("PARMS")=PARMS); 


  delete(mcmcAnal.Y); 
  delete(mcmcAnal.n_group); 
  delete(mcmcAnal.sd); 
  delete(mcmcAnal.doses); 
  delete(mcmcAnal.prior); 
  delete(output.parms); 
  delete(output.BMDS); 

  return wrap(data_out);
}


 
