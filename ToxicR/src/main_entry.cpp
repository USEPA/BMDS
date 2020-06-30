#include <RcppEigen.h>
#include <RcppGSL.h>

#include <iostream>
#include <string>
#include <vector>

#include "bmds_entry.h"


#include  <statmod.h>

#include <log_likelihoods.h>
#include <normal_likelihoods.h>
#include <normalModels.h>
#include <binomModels.h>
#include <IDPrior.h>

#include "bmd_calculate.h"

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
#include "dichotomous_entry_code.h"
#include "list_r_conversion.h"

using namespace Rcpp;
using namespace std;

#define MAX_PARMS 32 // Should never get close to this many!!!

// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::depends(RcppEigen)]]
//////////////////////////////////////////////////////////////////////////
// function: run_single_dichotomous
// purpose: takes input, which is assumed to be correct (i.e., filtered
// correctly by the R calling function), and then calls the library to
// run the corresponding analysis.
// output: BMD analysis with the model specified by NumericVector model
//
// [[Rcpp::export]]
List run_single_dichotomous(NumericVector model,
                            Eigen::MatrixXd data, Eigen::MatrixXd pr,
                            NumericVector options1, IntegerVector options2) 
{

  dichotomous_analysis Anal; 
  Anal.BMD_type =  (options1[0]==1)?eExtraRisk:eAddedRisk;
  Anal.BMR      =  options1[0]; 
  Anal.alpha    = options1[1];
  Anal.parms    = pr.rows(); 
  Anal.model    = (dich_model)model[0]; 
  Anal.Y        = new double[data.rows()] ; 
  Anal.n_group  = new double[data.rows()] ; 
  Anal.doses    = new double[data.rows()] ; 
  Anal.prior    = new double[pr.cols()*pr.rows()];
  Anal.prior_cols = pr.cols(); 
  Anal.n          = data.rows(); 
  Anal.degree = options2[1]; 

  if (Anal.model == dich_model::d_multistage){
    Anal.degree = Anal.parms - 1; 
  }

  for (int i = 0; i < data.rows(); i++){
    Anal.Y[i] = data(i,1); 
    Anal.n_group[i] = data(i,2); 
  }
  
  for (int i = 0; i < data.rows(); i++){
    Anal.doses[i] = data(i,0); 
  }
  
  // copy in column major order  
  for (int i = 0; i < pr.rows(); i++){
    for (int j = 0; j < pr.cols(); j++){
      Anal.prior[i + j*pr.rows()] = pr(i,j); 
    }
  }
  
  dichotomous_model_result res; 
  res.parms = new double[pr.rows()]; 
  res.cov   = new double[pr.rows()*pr.rows()]; 
  res.dist_numE = 200; 
  res.bmd_dist = new double[res.dist_numE*2]; 
  
  estimate_sm_laplace(&Anal, &res); 
  List rV = convert_dichotomous_fit_to_list(&res); 
  
  delete(Anal.Y); 
  delete(Anal.n_group); 
  delete(Anal.doses); 
  delete(Anal.prior); 
  delete(res.parms);   
  delete(res.cov);    
  delete(res.bmd_dist);  
  return rV;
}
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
// function: List run_ma_dichotomous()
// Purpose:  runs a model average based on the prior
//
// [[Rcpp::export]]
List run_ma_dichotomous(Eigen::MatrixXd data, List priors, NumericVector model_p,
                        NumericVector options1, IntegerVector options2){
	gsl_set_error_handler_off();
		///////////////////////////////////////////////////////
		//set up the MA_priors to send to the analysis
		MA_PRIORS ma_priors;
		ma_priors.mean_logistic = new(double[2]); ma_priors.sd_logistic = new(double[2]);
		ma_priors.mean_probit   = new(double[2]); ma_priors.sd_probit   = new(double[2]);
		ma_priors.mean_loglogit = new(double[3]); ma_priors.sd_loglogit = new(double[3]);
		ma_priors.mean_logprobit= new(double[3]); ma_priors.sd_logprobit= new(double[3]);
		ma_priors.mean_weibull  = new(double[3]); ma_priors.sd_weibull  = new(double[3]);
		ma_priors.mean_gamma    = new(double[3]); ma_priors.sd_gamma    = new(double[3]);
		ma_priors.mean_mult2    = new(double[3]); ma_priors.sd_mult2    = new(double[3]);
		ma_priors.mean_qlinear  = new(double[3]); ma_priors.sd_qlinear  = new(double[3]);
		ma_priors.mean_hill     = new(double[4]); ma_priors.sd_hill     = new(double[4]);
		////////////////////////////////////////////////////////
		NumericMatrix md = wrap(priors[0]);
		ma_priors.mean_logistic[0] = md(0,1); ma_priors.mean_logistic[1] = md(1,1);
		ma_priors.sd_logistic[0]   = md(0,2); ma_priors.sd_logistic[1] = md(1,2);

		md = wrap(priors[1]);
		ma_priors.mean_probit[0] = md(0,1); ma_priors.mean_probit[1] = md(1,1);
		ma_priors.sd_probit[0]   = md(0,2); ma_priors.sd_probit[1] = md(1,2);

		md = wrap(priors[2]);
		ma_priors.mean_loglogit[0] = md(0,1); ma_priors.mean_loglogit[1] = md(1,1); ma_priors.mean_loglogit[2] = md(2,1);
		ma_priors.sd_loglogit[0]   = md(0,2); ma_priors.sd_loglogit[1]   = md(1,2);   ma_priors.sd_loglogit[2] = md(2,2);

		md = wrap(priors[3]);
		ma_priors.mean_logprobit[0] = md(0,1); ma_priors.mean_logprobit[1] = md(1,1); ma_priors.mean_logprobit[2] = md(2,1);
		ma_priors.sd_logprobit[0]   = md(0,2); ma_priors.sd_logprobit[1]   = md(1,2);   ma_priors.sd_logprobit[2] = md(2,2);

		md = wrap(priors[4]);
		ma_priors.mean_weibull[0] = md(0,1); ma_priors.mean_weibull[1] = md(1,1); ma_priors.mean_weibull[2] = md(2,1);
		ma_priors.sd_weibull[0]   = md(0,2); ma_priors.sd_weibull[1]   = md(1,2);   ma_priors.sd_weibull[2] = md(2,2);

		md = wrap(priors[5]);
		ma_priors.mean_gamma[0] = md(0,1); ma_priors.mean_gamma[1] = md(1,1); ma_priors.mean_gamma[2] = md(2,1);
		ma_priors.sd_gamma[0]   = md(0,2); ma_priors.sd_gamma[1]   = md(1,2);   ma_priors.sd_gamma[2] = md(2,2);

		md = wrap(priors[6]);
		ma_priors.mean_mult2[0] = md(0,1); ma_priors.mean_mult2[1] = md(1,1); ma_priors.mean_mult2[2] = md(2,1);
		ma_priors.sd_mult2[0]   = md(0,2); ma_priors.sd_mult2[1]   = md(1,2);   ma_priors.sd_mult2[2] = md(2,2);

		md = wrap(priors[7]);
		ma_priors.mean_qlinear[0] = md(0,1); ma_priors.mean_qlinear[1] = md(1,1);
		ma_priors.sd_qlinear[0]   = md(0,2); ma_priors.sd_qlinear[1]   = md(1,2);


		md = wrap(priors[8]);
		ma_priors.mean_hill[0] = md(0,1); ma_priors.mean_hill[1] = md(1,1); ma_priors.mean_hill[2] = md(2,1);
		ma_priors.mean_hill[3] = md(3,1);
		ma_priors.sd_hill[0]   = md(0,2); ma_priors.sd_hill[1]   = md(1,2);   ma_priors.sd_hill[2] = md(2,2);
		ma_priors.sd_hill[3]   = md(3,2);

		////////////////////////////////////////////////////////////////
  	vector<BMDSInputData_t> vInputData;
  	BMDS_D_Opts1_t opt_one; opt_one.bmr = options1[0]; opt_one.alpha = options1[1];
		opt_one.background = options1[2];
		BMDS_D_Opts2_t opt_two; opt_two.bmrType = (options2[0]==1)?eExtraRisk:eAddedRisk; // 1 if extra 2 if added.
		opt_two.degree = options2[1];

		//copy the data in
		int n = data.rows();
		BMDSInputData_t *zTemp = new BMDSInputData_t[n];

		for (int i = 0; i < n; i++){
		  zTemp[i].dose      = data(i,0);
		  zTemp[i].response  = data(i,1);
		  zTemp[i].groupSize = data(i,2);
		//  vInputData.push_back(zTemp);
		}
		////////////////////////////////////////////////////////////////
		// ANALYSIS

	    double mp[9]; double posterior[9];

	  for (int i =0 ; i<9; i++){ mp[i] = model_p[i];}
		double ma_bmd[3], bmd[9], bmdl[9], bmdu[9];
		bmd_MA( eDich_3 , zTemp, &ma_priors, mp,
				&opt_one, &opt_two,
				n, posterior,
				ma_bmd,bmd,bmdl,bmdu);


		////////////////////////////////////////////////////////////////
		// CLEANUP
		delete(zTemp);
		delete(ma_priors.mean_logistic); delete(ma_priors.sd_logistic);
		delete(ma_priors.mean_probit)  ; delete(ma_priors.sd_probit);
		delete(ma_priors.mean_loglogit); delete(ma_priors.sd_loglogit);
		delete(ma_priors.mean_logprobit);delete(ma_priors.sd_logprobit);
		delete(ma_priors.mean_weibull);  delete(ma_priors.sd_weibull);
		delete(ma_priors.mean_gamma);    delete(ma_priors.sd_gamma);
		delete(ma_priors.mean_mult2);    delete(ma_priors.sd_mult2);
		delete(ma_priors.mean_qlinear);  delete(ma_priors.sd_qlinear);
		delete(ma_priors.mean_hill);     delete(ma_priors.sd_hill);
		////////////////////////////////////////////////////////////////
		NumericVector  BMD     = NumericVector::create(ma_bmd[0],ma_bmd[1],ma_bmd[2]);
		NumericVector  POSTW   = NumericVector::create(posterior[0],posterior[1],
											           posterior[2],posterior[3],
											           posterior[4],posterior[5],
											           posterior[6],posterior[7],
											           posterior[8]);

  	    List data_out  = List::create( Named("BMD")= BMD,
		   							   Named("POSTERIORS") = POSTW);
									  // Named("F_INFO") = NumericVector::create(a.MAP,logP));  // the MAP ESTIMATE and the
																								 // Laplace Integrating Constant
		return data_out;

}





// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::depends(RcppEigen)]]
//////////////////////////////////////////////////////////////////////////
// function: run_single_continuous
// purpose: takes input, which is assumed to be correct (i.e., filtered
// correctly by the R calling function), and then calls the library to
// run the corresponding analysis.
// output: BMD analysis with the model specified by NumericVector model
// [[Rcpp::export]]
List run_continuous_single(IntegerVector model, 
                           Eigen::MatrixXd Y, Eigen::MatrixXd X,
                           Eigen::MatrixXd prior, NumericVector options,
                           IntegerVector dist_type){
    
    bool   is_increasing = (bool)options[4]; 	double alpha = (double)options[3];
    double tail_p = (double)options[2]; 	double bmrf  = (double)options[1];
    int    riskType = (int)options[0];   
    unsigned int samples = (unsigned int) options[5];
    
    ////////////////////////////////////////////////
    /// Set up the analysis
    ////////////////////////////////////////////////
    continuous_analysis anal; 
    anal.Y       =    new double[Y.rows()]; 
    anal.n       =    Y.rows(); 
    anal.n_group =    new double[Y.rows()]; 
    anal.sd      =    new double[Y.rows()]; 
    anal.doses   =    new double[Y.rows()]; 
    anal.model   =    (cont_model) model[0]; 
    anal.disttype     = dist_type[0]; 
    anal.isIncreasing = is_increasing; 
    anal.alpha        = 0.005; //alpha for analyses; 
    anal.BMD_type     = riskType; 
    anal.BMR          = bmrf; 
    anal.samples      = samples; 
    anal.tail_prob    = tail_p; 
    anal.suff_stat    = Y.cols()==3;
    anal.parms        = prior.rows();
    anal.prior_cols   = prior.cols(); 
    anal.prior   = new double[prior.rows()*prior.cols()]; 
    cp_prior(prior,anal.prior);
    
    for (int i = 0; i < Y.rows(); i++){
      anal.Y[i] = Y(i,0); 
      anal.doses[i] = X(i,0); 
      if (Y.cols() == 3){ //sufficient statistics
        anal.n_group[i] = Y(i,2);
        anal.sd[i]      = Y(i,1); 
      }
    }
    ////////////////////////////////////
    continuous_model_result *result = new_continuous_model_result( anal.model,
                                                                   anal.parms,
                                                                   200); //have 200 equally spaced values
    ////////////////////////////////////
    estimate_sm_laplace(&anal,result);
    List rV = convert_continuous_fit_to_list(result); 	
    // free up memory
    del_continuous_model_result(result); 
    del_continuous_analysis(anal);
    return rV; 
    
  }
