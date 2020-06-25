#ifndef _DICHOTOMOUS_ENTRY_CODE_H
#define _DICHOTOMOUS_ENTRY_CODE_H

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


#endif