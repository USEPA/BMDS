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


#ifndef _DICHOTOMOUS_ENTRY_CODE_H
#define _DICHOTOMOUS_ENTRY_CODE_H

void estimate_ma_MCMC(dichotomousMA_analysis *MA,
                      dichotomous_analysis   *DA,
                      dichotomousMA_result   *res,
                      ma_MCMCfits            *ma);

void estimate_ma_laplace(dichotomousMA_analysis *MA,
                         dichotomous_analysis *DA ,
                         dichotomousMA_result *res);

void estimate_sm_laplace(dichotomous_analysis *DA ,
                         dichotomous_model_result *res);

void estimate_sm_mcmc(dichotomous_analysis *DA, 
                      continuous_model_result *res,
                      bmd_analysis_MCMC *mcmc);

#endif