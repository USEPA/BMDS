#ifdef R_COMPILATION  
  #include <RcppEigen.h>
  #include <RcppGSL.h>
#else 
  #include <Eigen/Dense>
#endif

#include <gsl/gsl_randist.h>
#include "bmdStruct.h"

#include <string>
#include <vector>
#include <limits>
#include <math.h>
#include <stdio.h>
#include "bmds_entry.h"
#include "bmdStruct.h"
#include "continuous_model_functions.h"

using namespace Rcpp;
using namespace std;
using Eigen::Map;
using Eigen::MatrixXd;
using Rcpp::as;

//const Map<MatrixXd> A(as<Map<MatrixXd>>(AA));

#include  <statmod.h>

#include <log_likelihoods.h>
#include <normal_likelihoods.h>
#include <normalModels.h>
#include <binomModels.h>
#include <IDPrior.h>

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



#ifndef _CONTINUOUS_ENTRY_CODE_H
#define _CONTINUOUS_ENTRY_CODE_H

bool convertSStat(Eigen::MatrixXd Y, Eigen::MatrixXd X,
                  Eigen::MatrixXd *SSTAT, Eigen::MatrixXd *SSTAT_LN,
                  Eigen::MatrixXd *UX);
                  
void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove);
void removeCol(Eigen::MatrixXd& matrix, unsigned int colToRemove);

bmd_analysis laplace_logNormal(Eigen::MatrixXd Y,Eigen::MatrixXd X,
                               Eigen::MatrixXd prior, contbmd riskType, cont_model CM,
                               bool is_increasing, 
                               double bmrf,   double bk_prob, 
                               double alpha, double step_size,
                               Eigen::MatrixXd init = Eigen::MatrixXd::Zero(1,1));
                               
bmd_analysis laplace_Normal(Eigen::MatrixXd Y,Eigen::MatrixXd X,
                            Eigen::MatrixXd prior, contbmd riskType, cont_model CM,
                            bool is_increasing, bool bConstVar,
                            double bmrf,   double bk_prob, 
                            double alpha, double step_size,
                            Eigen::MatrixXd init = Eigen::MatrixXd::Zero(1,1));

void estimate_ma_MCMC(continuousMA_analysis *MA,
                      continuous_analysis *CA ,
                      continuousMA_result *res,
                      ma_MCMCfits         *ma);
                            
void transfer_continuous_model(bmd_analysis a, continuous_model_result *model); 

void bmd_range_find(continuousMA_result *res, 
					double *range);
					
void estimate_ma_laplace(continuousMA_analysis *MA,
                         continuous_analysis *CA ,
                         continuousMA_result *res);

void estimate_sm_laplace(continuous_analysis *CA ,
                         continuous_model_result *res);

void estimate_sm_mcmc(continuous_analysis *CA,
                      continuous_model_result *res,
                      bmd_analysis_MCMC *mcmc); 


#endif
