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
#include "mcmc_analysis.h"
#include "bmdStruct.h"


#ifndef _CONTINUOUS_ENTRY_CODE_H
#define _CONTINUOUS_ENTRY_CODE_H

bmd_analysis create_bmd_analysis_from_mcmc(unsigned int burnin, mcmcSamples s);
void transfer_mcmc_output(mcmcSamples a, bmd_analysis_MCMC *b); 

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

                            
void transfer_continuous_model(bmd_analysis a, continuous_model_result *model); 

void bmd_range_find(continuousMA_result *res, 
					double *range);
					

					
/* Function: estimate_ma_mcmc 
 * Purpose:  This function performs a continuous Model Average (MA) for dichotomous
 *           data.  This is done using MCMC estimation. 
 * Input  : 
 *  @continuousMA_analysis - Pointer to memory containing information for the MA
 *                             analysis
 *  @continuous_analysis   - Pointer to memory containing information about the basic
 *                             type of continuous analysis.  Here individual model informaiton
 *                             is discarded. 
 * Return Values:    
 *  @continuousMA_result   - Pointer to the MA result.  Inside the result are all of the 
 *                             individual fit information.  The fit information is returned
 *                             using a summary statistic format, and is computed using the
 *                             proper number of burnins. 
 *  @ma_MCMCfits            - Pointer to a structure that returns the individual MCMC results
 *   
 
 * SPECIAL NOTES:   1) All memory is assumed allocated by the calling function. 
 *                  As a result, any memory should be dealocated by the calling function. 
 *                  Not doing so will result in memory leaks, and may be problematic on 
 *                  server environments. 
 *                  2) These funcions are overloaded with continuous counterparts 
 *                  that use the exact same estimation algorithms for continuous data. 
 */					
void estimate_ma_MCMC(continuousMA_analysis    *MA ,
                           continuous_analysis *CA ,
                           continuousMA_result *res,
                           ma_MCMCfits         *ma);

/* Function: estimate_ma_mcmc 
 * Purpose:  This function performs a continuous Model Average (MA) for dichotomous
 *           data.  This is done using Laplace/MAP estimation. 
 * Input  : 
 *  @continuousMA_analysis - Pointer to memory containing information for the MA
 *                             analysis
 *  @continuous_analysis   - Pointer to memory containing information about the basic
 *                             type of continuous analysis.  Here individual model informaiton
 *                             is discarded. 
 * Return Values:    
 *  @continuousMA_result   - Pointer to the MA result.  Inside the result are all of the 
 *                             individual fit information.  The fit information is returned
 *                             using a summary statistic format, and is computed using the
 *                             proper number of burnins. 
 *                             
 * SPECIAL NOTES:   1) All memory is assumed allocated by the calling function. 
 *                  As a result, any memory should be dealocated by the calling function. 
 *                  Not doing so will result in memory leaks, and may be problematic on 
 *                  server environments. 
 *                  2) These funcions are overloaded with continuous counterparts 
 *                  that use the exact same estimation algorithms for continuous data. 
 */						
void estimate_ma_laplace(continuousMA_analysis *MA,
                         continuous_analysis *CA ,
                         continuousMA_result *res);

/* Function: estimate_sm_laplace
 * Purpose:  This function performs a single model estimate for continuous data  
 * This is done using Laplace estimation. 
 * Input  : 
 *
 *  @continous_analysis   - Pointer to memory containing information about the basic
 *                             type of continuous analysis. Unlike the MA one individual 
 *                             model informaiton needed for the analysis. 
 * Return Values:    
 *  @continuous_model_result   - Pointer to a single model. The fit information is returned
 *                             using a summary statistic format.  All statistics are computed using
 *                             methodologies described in Wheeler et al (2020).
 *  
 *   
 * SPECIAL NOTES:   1) All memory is assumed allocated by the calling function. 
 *                  As a result, any memory should be dealocated by the calling function. 
 *                  Not doing so will result in memory leaks! This may be very problematic
 *                  on server environments. 
 *                  2) These funcions are overloaded with continuous counterparts 
 *                  that use the exact same estimation algorithms for continuous data.  
 */
void estimate_sm_laplace(continuous_analysis *CA ,
                         continuous_model_result *res);

/* Function: estimate_sm_mcmc
 * Purpose:  This function performs a single model estimate for continuous data  
 * This is done using MCMC estimation. 
 * Input  : 
 *
 *  @continous_analysis   - Pointer to memory containing information about the basic
 *                             type of continuous analysis. Unlike the MA one individual 
 *                             model informaiton needed for the analysis. 
 * Return Values:    
 *  @continuous_model_result   - Pointer to a single model. The fit information is returned
 *                             using a summary statistic format.  All statistics are computed using
 *                             methodologies described in Wheeler et al (2020).
 *  
 *   
 * SPECIAL NOTES:   1) All memory is assumed allocated by the calling function. 
 *                  As a result, any memory should be dealocated by the calling function. 
 *                  Not doing so will result in memory leaks! This may be very problematic
 *                  on server environments. 
 *                  2) These funcions are overloaded with continuous counterparts 
 *                  that use the exact same estimation algorithms for continuous data.  
 */
void estimate_sm_mcmc(continuous_analysis *CA,
                      continuous_model_result *res,
                      bmd_analysis_MCMC  *mcmc); 


#endif
