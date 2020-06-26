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
                            Eigen::MatrixXd data, Eigen::MatrixXd prior,
                            NumericVector options1, IntegerVector options2) {

	gsl_set_error_handler_off();
    //Begin: Analysis Set up ///////////////////////////////////////////////
    DModelID_t ID;
    switch(int(model[0])){
    case 1:
            ID = eDHill;  break;
    case 2:
            ID = eGamma;  break;
    case 3:
            ID = eLogistic;break;
    case 4:
            ID = eLogLogistic;break;
    case 5:
            ID = eLogProbit;break;
    case 6:
            ID = eMultistage;break;
    case 7:
            ID = eProbit;break;
    case 8:
            ID = eQLinear;break;
    case 9:
            ID = eWeibull;break;
    default:
            ID = eWeibull; 
    }
    
    BMDSInputType_t type = eDich_3; // Specify a dichotomous run

    // set up the initial structure required for the analysis
    constexpr auto CDF_TABLE_SIZE = 99;
    DichotomousDeviance_t zDev;
    int n = data.rows();
    dGoF_t zGoF;  zGoF.pzRow = new GoFRow_t[n];
    BMD_ANAL returnV; returnV.PARMS = new(double[MAX_PARMS]);
    for (int i = 0; i<MAX_PARMS; i++){ returnV.PARMS[i] = 0;  }
    returnV.nparms = prior.rows(); 
    returnV.boundedParms = new bool[MAX_PARMS];
    returnV.deviance = &zDev;
    returnV.gof = &zGoF;
    returnV.aCDF = new double[CDF_TABLE_SIZE];
    returnV.nCDF = CDF_TABLE_SIZE;
    returnV.covM = new double[prior.rows() * prior.rows() ];

    vector<BMDSInputData_t> vInputData;

    BMDS_D_Opts1_t opt_one; opt_one.bmr = options1[0]; opt_one.alpha = options1[1];
    opt_one.background = options1[2];
    BMDS_D_Opts2_t opt_two; opt_two.bmrType = (options2[0]==1)?eExtraRisk:eAddedRisk; // 1 if extra 2 if added.
    opt_two.degree = options2[1];
    ////////////////////////////////////////////////////////////////////////////
    //copy the data in
    BMDSInputData_t zTemp;
    for (int i = 0; i < n; i++){
      zTemp.dose     = data(i,0);
      zTemp.response  = data(i,1);
      zTemp.groupSize = data(i,2);
      vInputData.push_back(zTemp);
    }
    ///////////////////////////////////////////////////////////////////////
    //copy the priors for the analysis
    PRIOR *prior_to_model = new PRIOR[prior.rows()];
    for (int i = 0; i < prior.rows(); i++){
      prior_to_model[i].type = prior(i,0);
      prior_to_model[i].initalValue = prior(i,1); 
      prior_to_model[i].stdDev = prior(i,2);
      prior_to_model[i].minValue = prior(i,3);
      prior_to_model[i].maxValue = prior(i,4);
    }
    //End: Anlysis Set up 
    ///////////////////////////////////////////////////////////////////////


    run_dmodel2(&ID,&returnV,&type,&vInputData[0],
                prior_to_model,&opt_one,&opt_two,&n);



    CharacterVector x = CharacterVector::create( "BMD", "BMDL" );
    ///////////////////////////////////////////////////////////////////////
    // general return info
    Eigen::VectorXd parms(returnV.nparms);
    for (int i = 0; i < returnV.nparms; i++){ parms[i] = returnV.PARMS[i]; };
    NumericVector  Options   = NumericVector::create( opt_one.bmr,opt_one.alpha,double(opt_two.bmrType),
                                                      double(opt_two.degree),opt_one.background);
    NumericVector  BMD       = NumericVector::create( returnV.BMD, returnV.BMDL ,returnV.BMDU);
    NumericVector  statistics = NumericVector::create(returnV.MAP,returnV.AIC);
   // setup deviance
	  Eigen::MatrixXd deviance_t(2,3);
	  deviance_t(0,0) = returnV.deviance->llFull; deviance_t(1,0) = returnV.deviance->llReduced;
	  deviance_t(0,1) = returnV.deviance->devFit; deviance_t(1,1) = returnV.deviance->devReduced;
	  deviance_t(0,2) = returnV.deviance->pvFit;  deviance_t(1,2) = returnV.deviance->pvReduced;
	  IntegerMatrix dev_df(3,2);
	  dev_df(0,0) = returnV.deviance->nparmFull; dev_df(1,0) = returnV.deviance->nparmFit;
	  dev_df(2,0) = returnV.deviance->nparmReduced;
	  dev_df(0,1) = 0; // zero degrees of freedom on complete model
	  dev_df(1,1) = returnV.deviance->dfFit; dev_df(2,1) = returnV.deviance->dfReduced;
	  List dev = List::create(Named("computed_deviance")=deviance_t,Named("df")=dev_df);
	  
	  //gof table
	  Eigen::MatrixXd gof_table(n,6);
	  for (int i = 0; i < n; i++){
	    gof_table(i,0) = returnV.gof->pzRow[i].dose;
	    gof_table(i,1) = returnV.gof->pzRow[i].estProb;
	    gof_table(i,2) = returnV.gof->pzRow[i].expected;
	    gof_table(i,3) = returnV.gof->pzRow[i].observed;
	    gof_table(i,4) = returnV.gof->pzRow[i].scaledResidual;
	    gof_table(i,5) = returnV.gof->pzRow[i].size;
	  }
	  
	  Eigen::MatrixXd covMat(prior.rows(),prior.rows()); 
	  for (int i = 0 ; i < prior.rows(); i++){
	    for(int j = 0; j < prior.rows(); j++){
	      covMat(i,j) = returnV.covM[i + j*prior.rows()]; 
	    }
	  }
	  
	  Eigen::MatrixXd return_cdf(CDF_TABLE_SIZE,2); 
	  for (int ii = 0; ii < CDF_TABLE_SIZE; ii++){
	      return_cdf(ii,0) = returnV.aCDF[ii];
	      return_cdf(ii,1) = double(ii+1.0)/double(CDF_TABLE_SIZE+1.0); 
	  }
     	
     List data_out  = List::create( Named("options")=Options,
                                    Named("bmd") = BMD,
                                    Named("parameters") = parms,
                                    Named("cov") = covMat, 
                                    Named("fit_statistics") = statistics,
                               	    Named("deviance") = dev,
                                    Named("gof") = gof_table,
                                    Named("cdf") = return_cdf);
    
    
    delete returnV.boundedParms;
    delete zGoF.pzRow;
    delete returnV.PARMS;
    delete[] prior_to_model;
    delete[] returnV.aCDF;
    delete returnV.covM; 
    
    return data_out;
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
		BMDSInputData_t *zTemp = new(BMDSInputData_t[n]);

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
