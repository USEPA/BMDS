


#ifdef R_COMPILATION
#ifndef _LIST_R_CONVERSION_H
#define _LIST_R_CONVERSION_H
/////////////////////////////////////////////////////////////////////
List covert_MCMC_fit_to_list(bmd_analysis_MCMC *a);
/////////////////////////////////////////////////////////////////////
List covert_continuous_fit_to_list(continuous_model_result *result);
/////////////////////////////////////////////////////////////////////
List convert_continuous_maresults_to_list(continuousMA_result *result);
/////////////////////////////////////////////////////////////////////
List covert_dichotomous_fit_to_list(dichotomous_model_result *result);
#endif
#endif 