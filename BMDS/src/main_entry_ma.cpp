#include <RcppEigen.h>
#include <RcppGSL.h>

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
#include "continuous_entry_code.h"

/*
 * 
 * 
 */
List covert_continuous_fit_to_list(continuous_model_result *result){
  NumericVector  parms(result->nparms); 
  NumericMatrix  covM(result->nparms,result->nparms); 
  
  for (int i = 0; i < result->nparms; i++){
	  parms[i] = result->parms[i]; 
	  for (int j = 0; j < result->nparms; j++){
		covM(i,j) = result->cov[i + j*result->nparms]; 
	  }
  } 
  char dist[80];
  char str[160]; 
  
  switch(result->dist){
	   case distribution::normal:
			sprintf(dist,"Distribution: %s","Normal");
	   break; 
	   case distribution::normal_ncv:
			sprintf(dist,"Distribution: %s","Normal-NCV");
	   break; 
	   case distribution::log_normal:
			sprintf(dist,"Distribution: %s","Log-Normal");
	   break; 
  }
  
  switch(result->model){
	case cont_model::hill:
		sprintf(str,"Model: %s %s", "Hill",dist); 
	break; 
	case cont_model::exp_3:
		sprintf(str,"Model: %s %s", "Exponential-3",dist); 
	break;
	case cont_model::exp_5:
		sprintf(str,"Model: %s %s", "Exponential-5",dist); 
	break;	
	case cont_model::power: 
		sprintf(str,"Model: %s %s", "Power",dist); 
	break;  
  }
  double maximum = result->max; 
  NumericMatrix bmd_distribution(result->dist_numE , 2);
  
  for (int i = 0; i < result->dist_numE; i++){
	bmd_distribution(i,0) = result->bmd_dist[i]; 
	bmd_distribution(i,1) = result->bmd_dist[i+result->dist_numE];  

  } 

  List rV = List::create(Named("full_model") = str,
						 Named("parameters") = parms, 
                         Named("covariance") = covM, 
                         Named("bmd_dist")   = bmd_distribution,
                         Named("maximum")    = maximum);  
  return rV; 

}
 
/////////////////////////////////////////////////////////////////////////////
//
//
//
//
////////////////////////////////////////////////////////////////////////////
List convert_continuous_maresults_to_list(continuousMA_result *result){
	
	List fittedModels; 
	char str[80]; 
	
	for(int i = 0; i < result->nmodels ; i++){
		sprintf(str,"Fitted_Model_%d",i+1);
		fittedModels.push_back(covert_continuous_fit_to_list(result->models[i]),
							   str); 
	}
	NumericMatrix ma_bmd_dist(result->dist_numE,2); 
	NumericVector post_probs(result->nmodels);
	for(int i = 0; i < result->dist_numE; i++){
		ma_bmd_dist(i,0) = result->bmd_dist[i]; 
		ma_bmd_dist(i,1) = result->bmd_dist[i + result->dist_numE];
	}
	for (int i = 0; i < result->nmodels; i++){
		post_probs[i] = result->post_probs[i]; 
	}
	
	fittedModels.push_back(ma_bmd_dist,"BMD_CDF"); 
	fittedModels.push_back(post_probs ,"posterior_probs"); 
	
	return fittedModels; 
}

/////////////////////////////////////////////////////////////////////////////
//
//
/////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List run_continuous_ma_laplace(List model_priors, NumericVector model_type, 
                               NumericVector dist_type,
      			                   Eigen::MatrixXd Y, Eigen::MatrixXd X,
      			                   NumericVector options){

     bool   is_increasing = (bool)options[4]; 	double alpha = (double)options[3];
	 double tail_p = (double)options[2]; 	double bmrf  = (double)options[1];
	 int    riskType = (int)options[0];   
     unsigned int samples = (unsigned int) options[5];
   
	 continuousMA_analysis ma_anal;
	
	 ma_anal.nmodels = model_priors.length(); 
	 ma_anal.modelPriors = new double [ma_anal.nmodels]; 
	 ma_anal.priors  = new double *[ma_anal.nmodels];
	 ma_anal.nparms  = new int [ma_anal.nmodels];
	 ma_anal.actual_parms = new int [ma_anal.nmodels];
	 ma_anal.prior_cols   = new int[ma_anal.nmodels];
	 ma_anal.models  = new int [ma_anal.nmodels];
	 ma_anal.disttype= new int [ma_anal.nmodels];
	 
	 continuousMA_result *ma_result = new continuousMA_result; 
	 ma_result->nmodels    = ma_anal.nmodels; 
	 ma_result->dist_numE  = 300; 
	 ma_result->bmd_dist   = new double[300*2]; 
	 ma_result->post_probs = new double[ma_anal.nmodels];
	 ma_result->models     = new continuous_model_result*[ma_anal.nmodels];
	 
	 
	 for (int i = 0; i < ma_anal.nmodels; i++){
		   ma_anal.modelPriors[i] = 1.0/double(ma_anal.nmodels); 
		   Eigen::MatrixXd temp = model_priors[i];
		   ma_anal.priors[i]    = new double[temp.rows()*temp.cols()];
		 
		   cp_prior( temp, ma_anal.priors[i]);
		   ma_anal.nparms[i]     = temp.rows(); 
		   ma_anal.prior_cols[i] = temp.cols();
		   ma_anal.models[i]     = (int) model_type[i]; 
		   ma_anal.disttype[i]   = (int) dist_type[i]; 
		   ma_result->models[i] = new_continuous_model_result( ma_anal.models[i],
															   ma_anal.nparms[i],
															   200); //have 200 equally spaced values
															   
	 }
   
	 /// Set up the other info
	 continuous_analysis anal; 
	 anal.Y       =    new double[Y.rows()]; 
	 anal.n       =    Y.rows(); 
	 anal.n_group =    new double[Y.rows()]; 
	 anal.sd      =    new double[Y.rows()]; 
	 anal.doses   =    new double[Y.rows()]; 
	 anal.prior   =    NULL;
	 anal.isIncreasing = is_increasing; 
	 anal.alpha        = 0.005; //alpha for analyses; 
	 anal.BMD_type     = riskType; 
	 anal.BMR          = bmrf; 
	 anal.samples      = samples; 
	 anal.tail_prob    = tail_p; 
	 anal.suff_stat    = Y.cols()==3;
	 
	 
	 for (int i = 0; i < Y.rows(); i++){
	   anal.Y[i] = Y(i,0); 
	   anal.doses[i] = X(i,0); 
	   if (Y.cols() == 3){ //sufficient statistics
	     anal.n_group[i] = Y(i,2);
	     anal.sd[i]      = Y(i,1); 
	   }
	 }

   estimate_ma_laplace(&ma_anal,&anal,ma_result);
   List rV = convert_continuous_maresults_to_list(ma_result); 	

   // free up memory
   for (unsigned int i = 0; i < ma_result->nmodels; i++){
	   del_continuous_model_result(ma_result->models[i]); 
   }
   
   delete ma_result->post_probs; 
   delete ma_result->bmd_dist; 
   delete ma_result; 
   del_continuous_analysis(anal);
   del_continuousMA_analysis(ma_anal); 
   return rV; 

}


List covert_MCMC_fit_to_list(bmd_analysis_MCMC *a){
  List rV; 
  NumericMatrix parameters(a->samples,a->nparms);  
  NumericMatrix BMDS(a->samples,1); 
  
  for (unsigned int i = 0; i < a->samples; i++){
    BMDS[i] = a->BMDS[i]; 
    for( unsigned int j = 0; j < a->nparms; j++){
       parameters(i,j) = a->parms[i +j*a->samples]; 
    }
  }
  rV = List::create(Named("BMD_samples")=BMDS,Named("PARM_samples")=parameters); 
  return rV; 
}


List convert_mcmc_results(const ma_MCMCfits *a){
  List rV; 
  char str[80]; 

  for (unsigned int i=0; i < a->nfits; i++){
    sprintf(str,"Fitted_Model_%d",i+1);
    rV.push_back(covert_MCMC_fit_to_list(a->analyses[i]),
                 str); 
    
  }
  return rV; 
}
/////////////////////////////////////////////////////////////////////////////
//
//
/////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List run_continuous_ma_mcmc(List model_priors, NumericVector model_type, 
                               NumericVector dist_type,
                               Eigen::MatrixXd Y, Eigen::MatrixXd X,
                               NumericVector options){
  unsigned int burnin = (unsigned int) options[6];
  bool   is_increasing = (bool)options[4]; 	double alpha = (double)options[3];
  double tail_p = (double)options[2]; 	double bmrf  = (double)options[1];
  int    riskType = (int)options[0];   
  unsigned int samples = (unsigned int) options[5];
  
  continuousMA_analysis ma_anal;
  
  ma_anal.nmodels = model_priors.length(); 
  ma_anal.modelPriors = new double [ma_anal.nmodels]; 
  ma_anal.priors  = new double *[ma_anal.nmodels];
  ma_anal.nparms  = new int [ma_anal.nmodels];
  ma_anal.actual_parms = new int [ma_anal.nmodels];
  ma_anal.prior_cols   = new int[ma_anal.nmodels];
  ma_anal.models  = new int [ma_anal.nmodels];
  ma_anal.disttype= new int [ma_anal.nmodels];
  
  continuousMA_result *ma_result = new continuousMA_result; 
  ma_MCMCfits model_mcmc_info; 
  model_mcmc_info.analyses = new bmd_analysis_MCMC*[ma_anal.nmodels]; 
  model_mcmc_info.nfits = ma_anal.nmodels; 
  ma_result->nmodels    = ma_anal.nmodels; 
  ma_result->dist_numE  = 300; 
  ma_result->bmd_dist   = new double[300*2]; 
  ma_result->post_probs = new double[ma_anal.nmodels];
  ma_result->models     = new continuous_model_result*[ma_anal.nmodels];
  
  
  for (int i = 0; i < ma_anal.nmodels; i++){
    ma_anal.modelPriors[i] = 1.0/double(ma_anal.nmodels); 
    Eigen::MatrixXd temp = model_priors[i];
    ma_anal.priors[i]    = new double[temp.rows()*temp.cols()];
    //cout << temp << endl; 
    cp_prior( temp, ma_anal.priors[i]);
    ma_anal.nparms[i]     = temp.rows(); 
    ma_anal.prior_cols[i] = temp.cols();
    ma_anal.models[i]     = (int) model_type[i]; 
    ma_anal.disttype[i]   = (int) dist_type[i]; 
    //cout << ma_anal.models[i] << " " << dist_type[i] << endl; 
    ma_result->models[i] = new_continuous_model_result( ma_anal.models[i],
                                                        ma_anal.nparms[i],
                                                                      200); //have 200 equally spaced values
    model_mcmc_info.analyses[i] = new_mcmc_analysis(ma_anal.models[i],
                                                    ma_anal.nparms[i],
                                                              samples);
    
    
  }
  
  /// Set up the other info
  continuous_analysis anal; 
  anal.Y       =    new double[Y.rows()]; 
  anal.n       =    Y.rows(); 
  anal.n_group =    new double[Y.rows()]; 
  anal.sd      =    new double[Y.rows()]; 
  anal.doses   =    new double[Y.rows()]; 
  anal.prior   =    NULL;
  anal.isIncreasing = is_increasing; 
  anal.alpha        = 0.005; //alpha for analyses; 
  anal.BMD_type     = riskType; 
  anal.BMR          = bmrf; 
  anal.samples      = samples;
  anal.burnin       = burnin; 
  anal.tail_prob    = tail_p; 
  anal.suff_stat    = Y.cols()==3;
  
  
  for (int i = 0; i < Y.rows(); i++){
    anal.Y[i] = Y(i,0); 
    anal.doses[i] = X(i,0); 
    if (Y.cols() == 3){ //sufficient statistics
      anal.n_group[i] = Y(i,2);
      anal.sd[i]      = Y(i,1); 
    }
  }
  
  estimate_ma_MCMC(&ma_anal,&anal,ma_result,&model_mcmc_info);
 
  List rV = convert_mcmc_results(&model_mcmc_info); 	
 // List t2 = convert_continuous_maresults_to_list(ma_result); 
 // List rv(Named("mcm
  //////////////////////////////////////////////////////////
  // free up memory
  for (unsigned int i = 0; i < ma_result->nmodels; i++){
    del_continuous_model_result(ma_result->models[i]); 
    del_mcmc_analysis(model_mcmc_info.analyses[i]);
  }
  
  delete ma_result->post_probs; 
  delete ma_result->bmd_dist; 
  delete ma_result; 
  del_continuous_analysis(anal);
  del_continuousMA_analysis(ma_anal); 
  return rV; 
  
}

