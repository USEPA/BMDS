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
#include "continuous_entry_code.h"

#ifndef _DICHOTOMOUS_ENTRY_CODE_H
#define _DICHOTOMOUS_ENTRY_CODE_H

/* Function: estimate_ma_mcmc 
 * Purpose:  This function performs a dichotomous Model Average (MA) for dichotomous
 *           data.  This is done using MCMC estimation. 
 * Input  : 
 *  @dichotomousMA_analysis - Pointer to memory containing information for the MA
 *                             analysis
 *  @dichotomous_analysis   - Pointer to memory containing information about the basic
 *                             type of dichotomous analysis.  Here individual model informaiton
 *                             is discarded. 
 * Return Values:    
 *  @dichotomousMA_result   - Pointer to the MA result.  Inside the result are all of the 
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
void estimate_ma_MCMC(dichotomousMA_analysis *MA,
                      dichotomous_analysis   *DA,
                      dichotomousMA_result   *res,
                      ma_MCMCfits            *ma);

/* Function: estimate_ma_mcmc 
 * Purpose:  This function performs a dichotomous Model Average (MA) for dichotomous
 *           data.  This is done using Laplace estimation. 
 * Input  : 
 *  @dichotomousMA_analysis - Pointer to memory containing information for the MA
 *                             analysis
 *  @dichotomous_analysis   - Pointer to memory containing information about the basic
 *                             type of dichotomous analysis.  Here individual model informaiton
 *                             is discarded. 
 * Return Values:    
 *  @dichotomousMA_result   - Pointer to the MA result.  Inside the result are all of the 
 *                             individual fit information.  The fit information is returned
 *                             using a summary statistic format.  This is computed using the 
 *                             methodologies described in Wheeler et al ()2020).
 *  
 *   
 * SPECIAL NOTES:   1) All memory is assumed allocated by the calling function. 
 *                  As a result, any memory should be dealocated by the calling function. 
 *                  Not doing so will result in memory leaks! This may be very problematic
 *                  on server environments. 
 *                  2) These funcions are overloaded with continuous counterparts 
 *                  that use the exact same estimation algorithms for continuous data.  
 */
void estimate_ma_laplace(dichotomousMA_analysis *MA,
                         dichotomous_analysis *DA ,
                         dichotomousMA_result *res);

/* Function: estimate_sm_laplace
 * Purpose:  This function performs a single model estimate for dichotomous data  This is done using Laplace estimation. 
 * Input  : 
 *
 *  @dichotomous_analysis   - Pointer to memory containing information about the basic
 *                             type of dichotomous analysis. Unlike the MA one individual 
 *                             model informaiton needed for the analysis. 
 * Return Values:    
 *  @dichotomous_model_result   - Pointer to a single model. The fit information is returned
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
void estimate_sm_laplace(dichotomous_analysis *DA ,
                         dichotomous_model_result *res,
                         bool do_a_rescale = true);



/* Function: estimate_sm_mcmc
 * Purpose:  This function performs a single model estimate for dichotomous data  
 * This is done using MCMC estimation. 
 * Input  : 
 *
 *  @dichotomous_analysis   - Pointer to memory containing information about the basic
 *                             type of dichotomous analysis. Unlike the MA one individual 
 *                             model informaiton needed for the analysis. 
 * Return Values:    
 *  @dichotomous_model_result   - Pointer to a single model. The fit information is returned
 *                             using a summary statistic format.  All statistics are computed using
 *                             methodologies described in Wheeler et al (2020).
 *  
 *  @bmd_analysis_MCMC          - All of the posterior samples from the MCMC analysis. 
 *  
 *   
 * SPECIAL NOTES:   1) All memory is assumed allocated by the calling function. 
 *                  As a result, any memory should be dealocated by the calling function. 
 *                  Not doing so will result in memory leaks! This may be very problematic
 *                  on server environments. 
 *                  2) These funcions are overloaded with continuous counterparts 
 *                  that use the exact same estimation algorithms for continuous data.  
 */
void estimate_sm_mcmc(dichotomous_analysis *DA, 
                      dichotomous_model_result *res,
                      bmd_analysis_MCMC *mcmc,
                      bool do_a_rescale = true);

#endif