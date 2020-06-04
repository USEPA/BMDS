/*
 * continuous_entry_code.cpp
 * 
 * Copyright 2020  NIEHS <matt.wheeler@nih.gov>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */


 /* 
// a column order matrix px5 where p is 
// the number of parametersd
*/

#include "continuous_entry_code.h"
#include "mcmc_analysis.h"
#include <algorithm>
#include <vector>
#include <limits>

using namespace std; 
bool convertSStat(Eigen::MatrixXd Y, Eigen::MatrixXd X,
                  Eigen::MatrixXd *SSTAT, Eigen::MatrixXd *SSTAT_LN,
                  Eigen::MatrixXd *UX){
  bool convert = true; 
  
  if (Y.cols() == 1 ){
    // check to see if it can be converted into sufficient statistics4
    int temp = 0; 
    for (int i = 0; i < X.rows(); i++){
      for (int j = 0 ; j < X.rows(); j++){
        if (X(i,0) == X(j,0)){
          temp++; 
        }
      }
      if( temp == 1){
        convert = false; 
      }
    }
    
    if (convert){
      // we can convert the data
       *SSTAT    = createSuffStat( Y, X, false);
       *SSTAT_LN = createSuffStat(Y,X,true); 
        std::vector<double> uniqueX = unique_list(X );
        Eigen::MatrixXd temp_X(uniqueX.size(),1);
        for (int i = 0; i < uniqueX.size(); i++){
          temp_X(i,0) = uniqueX[i]; 
        }
        *UX = temp_X; 
        return true; 
 
     }
  
  }
  
  
    
  return false; 
}
void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove)
{
  unsigned int numRows = matrix.rows()-1;
  unsigned int numCols = matrix.cols();
  
  if( rowToRemove < numRows )
    matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);
  
  matrix.conservativeResize(numRows,numCols);
}

void removeCol(Eigen::MatrixXd& matrix, unsigned int colToRemove)
{
  unsigned int numRows = matrix.rows();
  unsigned int numCols = matrix.cols()-1;
  
  if( colToRemove < numCols )
    matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);
  
  matrix.conservativeResize(numRows,numCols);
}


bmd_analysis laplace_logNormal(Eigen::MatrixXd Y,Eigen::MatrixXd X,
                               Eigen::MatrixXd prior, contbmd riskType, cont_model CM,
                               bool is_increasing, 
                               double bmrf,   double bk_prob, 
                               double alpha, double step_size) {
  
  bool suff_stat = Y.cols() == 1? false:true; 
  
  std::vector<bool> fixedB(prior.rows());
  std::vector<double> fixedV(prior.rows());
  for (int i = 0; i < prior.rows(); i++) {
    fixedB[i] = false;
    fixedV[i] = 0.0;
  }
  
  
  IDcontinuousPrior model_prior(prior);
  
  lognormalEXPONENTIAL_BMD_NC likelihood_lnexp5U(Y, X, suff_stat,  NORMAL_EXP5_UP);
  lognormalEXPONENTIAL_BMD_NC likelihood_lnexp3U(Y, X, suff_stat,  NORMAL_EXP3_UP);
  lognormalEXPONENTIAL_BMD_NC likelihood_lnexp5D(Y, X, suff_stat,  NORMAL_EXP5_DOWN);
  lognormalEXPONENTIAL_BMD_NC likelihood_lnexp3D(Y, X, suff_stat,  NORMAL_EXP3_DOWN);
  lognormalHILL_BMD_NC likelihood_lnhill(Y, X, suff_stat,  0);
  bmd_analysis a;
  switch (CM)
  {
  case cont_model::hill:
    cout << "Running Exponential 4 Model Log-Normality Assumption." << endl;
    a = bmd_analysis_CNC<lognormalHILL_BMD_NC, IDcontinuousPrior>
                            (likelihood_lnhill,  model_prior, fixedB, fixedV,
                              riskType, bmrf, bk_prob,
                              is_increasing, alpha, step_size);
    break; 
  case cont_model::exp_3:
    cout << "Running Exponential 3 Model Log-Normality Assumption." << endl;
    if (is_increasing){
      a =  bmd_analysis_CNC<lognormalEXPONENTIAL_BMD_NC, IDcontinuousPrior>
                                    (likelihood_lnexp3U,  model_prior, fixedB, fixedV,
                                     riskType, bmrf,bk_prob,
                                     is_increasing, alpha, step_size);
      
    }else{
      a = bmd_analysis_CNC<lognormalEXPONENTIAL_BMD_NC, IDcontinuousPrior>
                                    (likelihood_lnexp3D,  model_prior, fixedB, fixedV,
                                     riskType, bmrf,bk_prob,
                                     is_increasing, alpha, step_size);
    }
    removeRow(a.COV,2);
    removeCol(a.COV,2);
    break; 
  case cont_model::exp_5:
  default: 
    cout << "Running Exponential 5 Model Log-Normality Assumption." << endl;
    if (is_increasing){
      a = bmd_analysis_CNC<lognormalEXPONENTIAL_BMD_NC, IDcontinuousPrior>
                                 (likelihood_lnexp5U,  model_prior, fixedB, fixedV,
                                  riskType, bmrf,bk_prob,
                                  is_increasing, alpha, step_size);
    }else{
      a = bmd_analysis_CNC<lognormalEXPONENTIAL_BMD_NC, IDcontinuousPrior>
                                  (likelihood_lnexp5D,  model_prior, fixedB, fixedV,
                                  riskType, bmrf,bk_prob,
                                  is_increasing, alpha, step_size);
    }
  break; 
  
  }
  
  return a; 
}


bmd_analysis laplace_Normal(Eigen::MatrixXd Y,Eigen::MatrixXd X,
                            Eigen::MatrixXd prior, contbmd riskType, cont_model CM,
                            bool is_increasing, bool bConstVar,
                            double bmrf,   double bk_prob, 
                            double alpha, double step_size) {
  
  bool suff_stat = Y.cols() == 1? false:true; 
  
  std::vector<bool> fixedB(prior.rows());
  std::vector<double> fixedV(prior.rows());
  
  for (int i = 0; i < prior.rows(); i++) {
    fixedB[i] = false;
    fixedV[i] = 0.0;
  }
  
  
  IDcontinuousPrior model_prior(prior);
  normalHILL_BMD_NC  likelihood_nhill(Y, X, suff_stat, bConstVar, 0);
  normalPOWER_BMD_NC likelihood_power(Y, X, suff_stat, bConstVar, 0);
  normalEXPONENTIAL_BMD_NC likelihood_nexp5U(Y, X, suff_stat, bConstVar, NORMAL_EXP5_UP);
  normalEXPONENTIAL_BMD_NC likelihood_nexp3U(Y, X, suff_stat, bConstVar, NORMAL_EXP3_UP);
  normalEXPONENTIAL_BMD_NC likelihood_nexp5D(Y, X, suff_stat, bConstVar, NORMAL_EXP5_DOWN);
  normalEXPONENTIAL_BMD_NC likelihood_nexp3D(Y, X, suff_stat, bConstVar, NORMAL_EXP3_DOWN);

  bmd_analysis a;
  switch (CM)
  {
  case cont_model::power:
    if (bConstVar){
      cout << "Running Power Model Normality Assumption using Laplace." << endl;
    }else{
      cout << "Running Power Model Normality-NCV Assumption using Laplace." << endl;
    } 
    a = bmd_analysis_CNC<normalPOWER_BMD_NC, IDcontinuousPrior>
                          (likelihood_power,  model_prior, fixedB, fixedV,
                           riskType, bmrf, bk_prob,
                           is_increasing, alpha, step_size);
    break; 
  case cont_model::hill:
    if (bConstVar){
      cout << "Running Hill Model Normality Assumption using Laplace." << endl;
    }else{
      cout << "Running Hill Model Normality-NCV Assumption using Laplace." << endl;
    }
    a = bmd_analysis_CNC<normalHILL_BMD_NC, IDcontinuousPrior>
                    (likelihood_nhill,  model_prior, fixedB, fixedV,
                     riskType, bmrf, bk_prob,
                     is_increasing, alpha, step_size);
    break; 
  case cont_model::exp_3:
    if (bConstVar){
      cout << "Running Exponential 3 Model Normality Assumption using Laplace." << endl;
    }else{
      cout << "Running Exponential 3 Model Normality-NCV Assumption using Laplace." << endl;
    }
    if (is_increasing){
      a =  bmd_analysis_CNC<normalEXPONENTIAL_BMD_NC, IDcontinuousPrior>
                          (likelihood_nexp3U,  model_prior, fixedB, fixedV,
                           riskType, bmrf,bk_prob,
                           is_increasing, alpha, step_size);
    }else{
      a = bmd_analysis_CNC<normalEXPONENTIAL_BMD_NC, IDcontinuousPrior>
                            (likelihood_nexp3D,  model_prior, fixedB, fixedV,
                             riskType, bmrf,bk_prob,
                             is_increasing, alpha, step_size);
    }
    removeRow(a.COV,2);
    removeCol(a.COV,2);
    break; 
  case cont_model::exp_5:
  default: 
    if (bConstVar){
      cout << "Running Exponential 5 Model Normality Assumption using Laplace." << endl;
    }else{
      cout << "Running Exponential 5 Model Normality-NCV Assumption using Laplace." << endl;
    }
    
    if (is_increasing){
      a = bmd_analysis_CNC<normalEXPONENTIAL_BMD_NC, IDcontinuousPrior>
                                    (likelihood_nexp5U,  model_prior, fixedB, fixedV,
                                     riskType, bmrf,bk_prob,
                                     is_increasing, alpha, step_size);
    }else{
      a = bmd_analysis_CNC<normalEXPONENTIAL_BMD_NC, IDcontinuousPrior>
                                    (likelihood_nexp5D,  model_prior, fixedB, fixedV,
                                    riskType, bmrf,bk_prob,
                                    is_increasing, alpha, step_size);
    }
    break; 
    
  }
  
  return a; 
}

void transfer_continuous_model(bmd_analysis a, continuous_model_result *model){
	if (model){
		model->max = a.MAP; 
		for (int i = 0; i< model->dist_numE; i ++){
				double temp = double(i)/double(model->dist_numE); 
				model->bmd_dist[i] = a.BMD_CDF.inv(temp);     // BMD @ probability
				model->bmd_dist[model->dist_numE + i] = temp; // probability 
		}
		for (int i = 0; i < model->nparms; i++){
				model->parms[i] = a.MAP_ESTIMATE(i,0); 
				for (int j = 0; j < model->nparms; j++){
					model->cov[i + j*model->nparms] = a.COV(i,j); 
				}
		}
	}
}


void bmd_range_find(continuousMA_result *res, 
					double *range){
 // assume the minimum BMD for the MA is always 0
 range[0] = 0.0; 
 double current_max = 0.0; 
 for (int j = 10; j > 1; j--){
	 for (int i = 0; i < res->nmodels;i++){
		int temp_idx = res->models[i]->dist_numE -j; 
		
		// make sure we are not dealing with an infinite value
		// or not a number
		if (!isnan(res->models[i]->bmd_dist[temp_idx]) && 
			!isinf(res->models[i]->bmd_dist[temp_idx])){
			if ( res->models[i]->bmd_dist[temp_idx] > current_max){
				current_max = res->models[i]->bmd_dist[temp_idx]; 
			}
		}
		 
	 }
 }
 // if we don't find a max then the upper limit is NAN
 range[1] = current_max == 0.0 ? std::numeric_limits<double>::quiet_NaN():current_max; 
  
}

void estimate_ma_laplace(continuousMA_analysis *MA,
                         continuous_analysis *CA ,
                         continuousMA_result *res){
  // standardize the data
  int n_rows = CA->n; int n_cols = CA->suff_stat?3:1; 
  //cerr << "Sufficient Stat: " << n_cols << endl; 
  Eigen::MatrixXd Y(n_rows,n_cols); 
  Eigen::MatrixXd X(n_rows,1); 
  // copy the origional data
  for (int i = 0; i < n_rows; i++){
    Y(i,0) = CA->Y[i]; 
    X(i,0) = CA->doses[i]; 
    if(CA->suff_stat){
      
      Y(i,1) = CA->sd[i]; 
      Y(i,2) = CA->n_group[i]; 
    }
  }

  Eigen::MatrixXd SSTAT, SSTAT_LN, UX; 
  Eigen::MatrixXd Y_LN, Y_N;
  if(!CA->suff_stat){
    //convert to sufficient statistics for speed if we can
     CA->suff_stat = convertSStat(Y, X, &SSTAT, &SSTAT_LN,&UX); 
    if (CA->suff_stat)// it can be converted
    {
      X = UX; 
      Y_N = cleanSuffStat(SSTAT,UX,false);  
      Y_LN = cleanSuffStat(SSTAT_LN,UX,true); 
      
    }
  }else{
      SSTAT = cleanSuffStat(Y,X,false); 
      SSTAT_LN = cleanSuffStat(Y,X,true);
      std::vector<double> tux = unique_list(X); 
      UX = Eigen::MatrixXd(tux.size(),1); 
      for (unsigned int i = 0; i < tux.size(); i++){
         UX(i,0) = tux[i]; 
      }
      Y_N = SSTAT; 
      X = UX; 
      Y_LN = SSTAT_LN; 
  }

  double divisor = get_diviosor( Y,  X); 
  double  max_dose = X.maxCoeff(); 
   
  if (CA->suff_stat){
    X = UX; 
  //  Y_N = cleanSuffStat(SSTAT,UX,false);  
  //  Y_LN = cleanSuffStat(SSTAT_LN,UX,true); 
     Eigen::MatrixXd temp; 
     temp = Y_N.col(2);
     Y_N.col(2) = Y_N.col(1);
     Y_N.col(1) = temp; 
     temp = Y_LN.col(2);
     Y_LN.col(2) = Y_LN.col(1);
     Y_LN.col(1) = temp; 
     X = X/max_dose;
  }else{
    Y = (1/divisor)*Y; // scale the data with the divisor term. 
    X = X/max_dose;
  }

  
  bmd_analysis a[MA->nmodels];
  
  for (int i = 0; i < MA->nmodels; i++ ){
   
    Eigen::MatrixXd tprior(MA->nparms[i],MA->prior_cols[i]);
    for (int m = 0; m < MA->nparms[i]; m++){
      for (int n = 0; n < MA->prior_cols[i]; n++){
        tprior(m,n) = MA->priors[i][m + n*MA->nparms[i]]; 
      }
    }
      
    if (MA->disttype[i] == distribution::log_normal){
      
      if (CA->suff_stat ){
        a[i] = laplace_logNormal(Y_LN, X,
                                  tprior, CA->BMD_type, (cont_model)MA->models[i],
                                  CA->isIncreasing, CA->BMR, 
                                  CA->tail_prob,  
                                  CA->alpha, 0.02);
      }else{
        a[i] = laplace_logNormal(Y, X,
                                 tprior, CA->BMD_type, (cont_model)MA->models[i],
                                 CA->isIncreasing, CA->BMR, 
                                 CA->tail_prob,  
                                 CA->alpha, 0.02);
        
      }
    
    }else{
      
      bool isNCV = MA->disttype[i] == distribution::normal_ncv? false:true; 
      if (CA->suff_stat ){
        a[i] = laplace_Normal(Y_N, X,
                              tprior, CA->BMD_type, (cont_model)MA->models[i],
                              CA->isIncreasing,isNCV, CA->BMR, 
                              CA->tail_prob,  
                              CA->alpha, 0.02);
      }else{
        a[i] = laplace_Normal(Y, X,
                              tprior, CA->BMD_type, (cont_model)MA->models[i],
                              CA->isIncreasing,isNCV, CA->BMR, 
                              CA->tail_prob,  
                              CA->alpha, 0.02);
        
      }
     
    }
//FIX ME: RESCALE THE FITS TO WHAT THEY SHOULD BE
  }
  
  
  double post_probs[MA->nmodels]; 
  double temp =0.0; 
  double max_prob = -1.0*std::numeric_limits<double>::infinity(); 
  for (int i = 0; i < MA->nmodels; i++){
        temp  = 	a[i].MAP_ESTIMATE.rows()/2 * log(2 * M_PI) - a[i].MAP + 0.5*log(max(0.0,a[i].COV.determinant()));
        max_prob = temp > max_prob? temp:max_prob; 
        post_probs[i] = temp; 
  }
  
  
  double norm_sum = 0.0; 
  
  
  for (int i = 0; i < MA->nmodels; i++){
        post_probs[i] = post_probs[i] - max_prob + log(MA->modelPriors[i]); //FIX ME: ADD MODEL PROBS
        norm_sum += exp(post_probs[i]);
        post_probs[i] = exp(post_probs[i]);
  }

  for (int j = 0; j < MA->nmodels; j++){
    post_probs[j] = post_probs[j]/ norm_sum; 
    for (double  i = 0.0; i < 0.99; i += 0.01 ){
      if ( isnan(a[j].BMD_CDF.inv(i))){
        post_probs[j] = 0; // if the cdf has nan in it then it needs a 0 posterior
                           // probability
      }  
      
    } 
  }
  norm_sum = 0.0; 
  for (int i =0; i < MA->nmodels; i++){
    norm_sum += post_probs[i]; 
  }
 
  for (int i =0; i < MA->nmodels; i++){
    post_probs[i] = post_probs[i]/norm_sum; 
    res->post_probs[i] = post_probs[i];
	
	transfer_continuous_model(a[i],res->models[i]);
	res->models[i]->model = MA->models[i]; 
	res->models[i]->dist  = MA->disttype[i];
  }
	 
  double range[2]; 

  bmd_range_find(res,range);
  double range_bmd = range[1] - range[0]; 
  for (int i = 0; i < res->dist_numE; i ++){
	  double cbmd = double(i)/double(res->dist_numE)*range_bmd; 
	  double prob = 0.0; 
	  
	  for (int j = 0; j < MA->nmodels; j++){
			prob += isnan(a[j].BMD_CDF.P(cbmd))?0:a[j].BMD_CDF.P(cbmd)*post_probs[j]; 
      }
	  res->bmd_dist[i] = cbmd; 
	  res->bmd_dist[i+res->dist_numE]  = prob;
  }
	
  return; 
}


/*############################################################################
 * 
 * 
 * 
 * 
##############################################################################*/
mcmcSamples mcmc_logNormal(Eigen::MatrixXd Y,Eigen::MatrixXd X,
                            Eigen::MatrixXd prior, contbmd riskType, cont_model CM,
                            bool is_increasing, 
                            double bmrf,   double bk_prob, 
                            double alpha,  int samples, int burnin) {
   
  bool suff_stat = Y.cols() == 1? false:true; 
  std::vector<bool> fixedB(prior.rows());
  std::vector<double> fixedV(prior.rows());
  for (int i = 0; i < prior.rows(); i++) {
    fixedB[i] = false;
    fixedV[i] = 0.0;
  }

  mcmcSamples a;
  int adverseR; 
  switch (CM)
  {
  case cont_model::hill:
    cout << "Running Hill Model Log-Normality Assumption using MCMC." << endl;
    
     a =  MCMC_bmd_analysis_CONTINUOUS_LOGNORMAL<lognormalHILL_BMD_NC, IDcontinuousPrior>
                                    (Y,  X, prior, fixedB, fixedV, is_increasing,
                                     bk_prob,suff_stat,bmrf, riskType, alpha,
                                     samples,0); 
    break; 
  case cont_model::exp_3:
      adverseR = is_increasing?NORMAL_EXP3_UP: NORMAL_EXP3_DOWN; 
      cout << "Running Exponential 3 Model Log-Normality Assumption using MCMC." << endl;
      a =  MCMC_bmd_analysis_CONTINUOUS_LOGNORMAL<lognormalEXPONENTIAL_BMD_NC, IDcontinuousPrior>
                                              (Y,  X, prior, fixedB, fixedV, is_increasing,
                                              bk_prob,suff_stat,bmrf, riskType, alpha,
                                              samples,adverseR);
    break; 
  case cont_model::exp_5:
  default: 
      adverseR = is_increasing?NORMAL_EXP5_UP: NORMAL_EXP5_DOWN; 
      cout << "Running Exponential 5 Model Log-Normality Assumption using MCMC." << endl;
      a =  MCMC_bmd_analysis_CONTINUOUS_LOGNORMAL<lognormalEXPONENTIAL_BMD_NC, IDcontinuousPrior>
                                              (Y,  X, prior, fixedB, fixedV, is_increasing,
                                               bk_prob,suff_stat,bmrf, riskType, alpha,
                                               samples,adverseR);
    break; 
    
  }
  
  return a; 
}

mcmcSamples mcmc_Normal(Eigen::MatrixXd Y,Eigen::MatrixXd X,
                         Eigen::MatrixXd prior, contbmd riskType, cont_model CM,
                         bool is_increasing, bool bConstVar,
                         double bmrf,   double bk_prob, 
                         double alpha, int samples,
                         int burnin) {
  
  bool suff_stat = Y.cols() == 1? false:true; 
  std::vector<bool> fixedB(prior.rows());
  std::vector<double> fixedV(prior.rows());
  
  for (int i = 0; i < prior.rows(); i++) {
    fixedB[i] = false;
    fixedV[i] = 0.0;
  }
  
  mcmcSamples a;
  int adverseR; 
  switch (CM)
  {
  case cont_model::hill:
    if (bConstVar){
      cout << "Running Hill Model Normality Assumption using MCMC." << endl;
    }else{
      cout << "Running Hill Model Normality-NCV Assumption using MCMC." << endl;
    }
    
      a =  MCMC_bmd_analysis_CONTINUOUS_NORMAL<normalHILL_BMD_NC, IDcontinuousPrior>
                                                (Y,  X, prior, fixedB, fixedV, is_increasing,
                                                 bk_prob,suff_stat,bmrf, riskType,bConstVar, alpha,
                                                 samples,0); 
    break; 
  case cont_model::exp_3:
    adverseR = is_increasing?NORMAL_EXP3_UP: NORMAL_EXP3_DOWN; 
    if (bConstVar){
      cout << "Running Exponential 3 Model Normality Assumption using MCMC." << endl;
    }else{
      cout << "Running Exponential 3 Model Normality-NCV Assumption using MCMC." << endl;
    }
    a =  MCMC_bmd_analysis_CONTINUOUS_NORMAL<normalEXPONENTIAL_BMD_NC, IDcontinuousPrior>
                                                (Y,  X, prior, fixedB, fixedV, is_increasing,
                                                 bk_prob,suff_stat,bmrf, riskType,bConstVar, alpha,
                                                 samples,adverseR);
    break; 
  case cont_model::exp_5:
  
    adverseR = is_increasing?NORMAL_EXP5_UP: NORMAL_EXP5_DOWN; 
    if (bConstVar){
      cout << "Running Exponential 5 Model Normality Assumption using MCMC." << endl;
    }else{
      cout << "Running Exponential 5 Model Normality-NCV Assumption using MCMC." << endl;
    }
    a =  MCMC_bmd_analysis_CONTINUOUS_NORMAL<normalEXPONENTIAL_BMD_NC, IDcontinuousPrior>
                                                (Y,  X, prior, fixedB, fixedV, is_increasing,
                                                 bk_prob,suff_stat,bmrf, riskType,bConstVar, alpha,
                                                 samples,0);
    break; 
  case cont_model::power:
  default:  
    if (bConstVar){
      cout << "Running Power Model Normality Assumption using MCMC." << endl;
    }else{
      cout << "Running Powrer Model Normality-NCV Assumption using MCMC." << endl;
    }
    a =  MCMC_bmd_analysis_CONTINUOUS_NORMAL<normalPOWER_BMD_NC, IDcontinuousPrior>
                                              (Y,  X, prior, fixedB, fixedV, is_increasing,
                                              bk_prob,suff_stat,bmrf, riskType,bConstVar, alpha,
                                              samples,adverseR);
    
  
  break; 
    
  }
  //convert a stuff
  
  //
  return a; 
}

bmd_analysis create_bmd_analysis_from_mcmc(unsigned int burnin, mcmcSamples s){
  bmd_analysis rV;
  rV.MAP          = s.map; 
  rV.MAP_ESTIMATE = s.map_estimate; 
  rV.COV          = s.map_cov; 
  std::vector<double> v(s.BMD.cols()-burnin); 
  for (int i = burnin; i < s.BMD.cols();i++){
	v[i] = s.BMD(0,i);   
  
  }
	
  std::vector<double>  prob;
  std::vector<double> bmd_q;  
  std::sort(v.begin(), v.end());
  for (double k = 0.005; k <= 0.995; k += 0.005){
		prob.push_back(k); 
		prob.push_back(v[int(k*(v.size()-burnin))]);
  }
	
  return rV; 
}

void transfer_mcmc_output(mcmcSamples a, bmd_analysis_MCMC *b){
  if (b){
    b->samples = a.samples.cols(); 
    b->nparms  = a.samples.rows(); 

    for(unsigned int i= 0; i < a.BMD.cols();  i++){
      b->BMDS[i] = a.BMD(0,i); 
      for(unsigned int j = 0; j < a.samples.rows(); j++){
        b->parms[i + j*a.BMD.cols()] = a.samples(j,i); 
      }
    }
  }
}


void fixRescaleLogNormal(Eigen::MatrixXd Y,Eigen::MatrixXd X, cont_model CM,
                         Eigen::MatrixXd prior,double max_dose,double divisor,
                         bool is_increasing) {
  
  bool suff_stat = Y.cols() == 1? false:true; 
  std::vector<bool> fixedB(prior.rows());
  std::vector<double> fixedV(prior.rows());
  

  for (int i = 0; i < prior.rows(); i++) {
      fixedB[i] = false;
      fixedV[i] = 0.0;
  }

  int adverseR; 
  switch (CM)
  {
  case cont_model::hill:
    
    
    cBMDModel<lognormalHILL_BMD_NC, IDcontinuousPrior>  model(likelihood, prior, fixedB, fixedV, isIncreasing);
    // <lognormalHILL_BMD_NC, IDcontinuousPrior>
    break; 
  case cont_model::exp_3:
    cBMDModel<lognormalHILL_BMD_NC, IDcontinuousPrior>  model(likelihood, prior, fixedB, fixedV, isIncreasing);
    adverseR = is_increasing?NORMAL_EXP3_UP: NORMAL_EXP3_DOWN; 
    // <lognormalEXPONENTIAL_BMD_NC, IDcontinuousPrior>

    break; 
  case cont_model::exp_5:
  default: 
    cBMDModel<lognormalHILL_BMD_NC, IDcontinuousPrior>  model(likelihood, prior, fixedB, fixedV, isIncreasing);
    adverseR = is_increasing?NORMAL_EXP5_UP: NORMAL_EXP5_DOWN; 
    // <lognormalEXPONENTIAL_BMD_NC, IDcontinuousPrior>
     
    break; 
    
  }
  
  return; 
}


/*
void fixRescaleNormal(Eigen::MatrixXd Y,Eigen::MatrixXd X,
                        Eigen::MatrixXd prior, contbmd riskType, cont_model CM,
                        bool is_increasing, bool bConstVar,
                        double bmrf,   double bk_prob, 
                        double alpha, int samples,
                        int burnin) {
  
  bool suff_stat = Y.cols() == 1? false:true; 
  std::vector<bool> fixedB(prior.rows());
  std::vector<double> fixedV(prior.rows());
  
  for (int i = 0; i < prior.rows(); i++) {
    fixedB[i] = false;
    fixedV[i] = 0.0;
  }
  
  mcmcSamples a;
  int adverseR; 
  switch (CM)
  {
  case cont_model::hill:
    if (bConstVar){
      cout << "Running Hill Model Normality Assumption using MCMC." << endl;
    }else{
      cout << "Running Hill Model Normality-NCV Assumption using MCMC." << endl;
    }
    
    a =  MCMC_bmd_analysis_CONTINUOUS_NORMAL<normalHILL_BMD_NC, IDcontinuousPrior>
      (Y,  X, prior, fixedB, fixedV, is_increasing,
       bk_prob,suff_stat,bmrf, riskType,bConstVar, alpha,
       samples,0); 
    break; 
  case cont_model::exp_3:
    adverseR = is_increasing?NORMAL_EXP3_UP: NORMAL_EXP3_DOWN; 
    if (bConstVar){
      cout << "Running Exponential 3 Model Normality Assumption using MCMC." << endl;
    }else{
      cout << "Running Exponential 3 Model Normality-NCV Assumption using MCMC." << endl;
    }
    a =  MCMC_bmd_analysis_CONTINUOUS_NORMAL<normalEXPONENTIAL_BMD_NC, IDcontinuousPrior>
      (Y,  X, prior, fixedB, fixedV, is_increasing,
       bk_prob,suff_stat,bmrf, riskType,bConstVar, alpha,
       samples,adverseR);
    break; 
  case cont_model::exp_5:
    
    adverseR = is_increasing?NORMAL_EXP5_UP: NORMAL_EXP5_DOWN; 
    if (bConstVar){
      cout << "Running Exponential 5 Model Normality Assumption using MCMC." << endl;
    }else{
      cout << "Running Exponential 5 Model Normality-NCV Assumption using MCMC." << endl;
    }
    a =  MCMC_bmd_analysis_CONTINUOUS_NORMAL<normalEXPONENTIAL_BMD_NC, IDcontinuousPrior>
      (Y,  X, prior, fixedB, fixedV, is_increasing,
       bk_prob,suff_stat,bmrf, riskType,bConstVar, alpha,
       samples,0);
    break; 
  case cont_model::power:
  default:  
    if (bConstVar){
      cout << "Running Power Model Normality Assumption using MCMC." << endl;
    }else{
      cout << "Running Powrer Model Normality-NCV Assumption using MCMC." << endl;
    }
    a =  MCMC_bmd_analysis_CONTINUOUS_NORMAL<normalPOWER_BMD_NC, IDcontinuousPrior>
      (Y,  X, prior, fixedB, fixedV, is_increasing,
       bk_prob,suff_stat,bmrf, riskType,bConstVar, alpha,
       samples,adverseR);
    
    
    break; 
    
  }
  //convert a stuff
  
  //
  return a; 
}
*/
void estimate_ma_MCMC(continuousMA_analysis *MA,
                      continuous_analysis   *CA,
                      continuousMA_result   *res,
                      ma_MCMCfits           *ma){ 
  // standardize the data
  int n_rows = CA->n; int n_cols = CA->suff_stat?3:1; 
  //cerr << "Sufficient Stat: " << n_cols << endl; 
  Eigen::MatrixXd Y(n_rows,n_cols); 
  Eigen::MatrixXd X(n_rows,1); 
  // copy the origional data
  for (int i = 0; i < n_rows; i++){
    Y(i,0) = CA->Y[i]; 
    X(i,0) = CA->doses[i]; 
    if(CA->suff_stat){
      
      Y(i,1) = CA->sd[i]; 
      Y(i,2) = CA->n_group[i]; 
    }
  }
  
  Eigen::MatrixXd SSTAT, SSTAT_LN, UX; 
  Eigen::MatrixXd Y_LN, Y_N;
  if(!CA->suff_stat){
    //convert to sufficient statistics for speed if we can
    CA->suff_stat = convertSStat(Y, X, &SSTAT, &SSTAT_LN,&UX); 
    if (CA->suff_stat)// it can be converted
    {
      X = UX; 
      Y_N = cleanSuffStat(SSTAT,UX,false);  
      Y_LN = cleanSuffStat(SSTAT_LN,UX,true); 
      
    }
  }else{
    SSTAT = cleanSuffStat(Y,X,false); 
    SSTAT_LN = cleanSuffStat(Y,X,true);
    std::vector<double> tux = unique_list(X); 
    UX = Eigen::MatrixXd(tux.size(),1); 
    for (unsigned int i = 0; i < tux.size(); i++){
      UX(i,0) = tux[i]; 
    }
    Y_N = SSTAT; 
    X = UX; 
    Y_LN = SSTAT_LN; 
  }
  
  double divisor = get_diviosor( Y,  X); 
  double  max_dose = X.maxCoeff(); 
  
  if (CA->suff_stat){
    X = UX; 
    //  Y_N = cleanSuffStat(SSTAT,UX,false);  
    //  Y_LN = cleanSuffStat(SSTAT_LN,UX,true); 
    Eigen::MatrixXd temp; 
    temp = Y_N.col(2);
    Y_N.col(2) = Y_N.col(1);
    Y_N.col(1) = temp; 
    temp = Y_LN.col(2);
    Y_LN.col(2) = Y_LN.col(1);
    Y_LN.col(1) = temp; 
    X = X/max_dose;
  }else{
    Y = (1/divisor)*Y; // scale the data with the divisor term. 
    X = X/max_dose;
  }
  
  
  mcmcSamples a[MA->nmodels];
  unsigned int samples = CA->samples; 
  unsigned int burnin  = CA->burnin;  
  for (int i = 0; i < MA->nmodels; i++ ){
    
    Eigen::MatrixXd tprior(MA->nparms[i],MA->prior_cols[i]);
    for (int m = 0; m < MA->nparms[i]; m++){
        for (int n = 0; n < MA->prior_cols[i]; n++){
          tprior(m,n) = MA->priors[i][m + n*MA->nparms[i]]; 
        }
    }
    
    if (MA->disttype[i] == distribution::log_normal){
      
      if (CA->suff_stat ){
        a[i] = mcmc_logNormal(Y_LN, X,
                              tprior, CA->BMD_type, (cont_model)MA->models[i],
                                                      CA->isIncreasing, CA->BMR, 
                                                      CA->tail_prob,  
                                                      CA->alpha, samples, burnin);
      }else{
        a[i] = mcmc_logNormal(Y, X,
                                 tprior, CA->BMD_type, (cont_model)MA->models[i],
                                                      CA->isIncreasing, CA->BMR, 
                                                      CA->tail_prob,  
                                                      CA->alpha,samples, burnin);
        
      }
      
    }else{
      
      bool isNCV = MA->disttype[i] == distribution::normal_ncv? false:true; 
      if (CA->suff_stat ){
        a[i] = mcmc_Normal(Y_N, X,
                           tprior, CA->BMD_type, (cont_model)MA->models[i],
                           CA->isIncreasing,isNCV, CA->BMR, 
                           CA->tail_prob,  
                           CA->alpha, samples, burnin);
      }else{
        a[i] = mcmc_Normal(Y, X,
                          tprior, CA->BMD_type, (cont_model)MA->models[i],
                          CA->isIncreasing,isNCV, CA->BMR, 
                          CA->tail_prob,  
                          CA->alpha,samples, burnin);
        
      }
      
    }

     for (int m = 0; m < samples; m++){
          a[i].BMD(0,m) = a[i].BMD(0,m)*max_dose; 
          Eigen::MatrixXd temp = rescale_parms(a[i].samples.col(i), 
                                               (cont_model)MA->models[i] ,
                                                max_dose, divisor, MA->disttype[i] == distribution::log_normal); 
                                                
         
          a[i].samples.col(m) = temp.transpose(); 
         
     }
     a[i].map_estimate = rescale_parms(a[i].map_estimate, 
                                       (cont_model)MA->models[i] ,
                                        max_dose, divisor, MA->disttype[i] == distribution::log_normal); 
     a[i].map_cov      = rescale_cov_matrix(a[i].map_cov,
										                    a[i].map_estimate, 
                                       (cont_model)MA->models[i] ,
                                        max_dose, divisor, MA->disttype[i] == distribution::log_normal);
  }

  bmd_analysis b[MA->nmodels]; 
  
  for (int i = 0; i < MA->nmodels; i++){
    b[i] = create_bmd_analysis_from_mcmc(burnin,a[i]);
  }
  
  double post_probs[MA->nmodels]; 
  double temp =0.0; 
  double max_prob = -1e300; 
  for (int i = 0; i < MA->nmodels; i++){
    temp  = 	b[i].MAP_ESTIMATE.rows()/2 * log(2 * M_PI) - b[i].MAP + 0.5*log(max(0.0,b[i].COV.determinant()));
    max_prob = temp > max_prob? temp:max_prob; 
    post_probs[i] = temp; 
  }
  double norm_sum = 0.0; 
  
  for (int i = 0; i < MA->nmodels; i++){
    post_probs[i] = post_probs[i] - max_prob + log(MA->modelPriors[i]); //FIX ME: ADD MODEL PROBS
    norm_sum += exp(post_probs[i]);
    post_probs[i] = exp(post_probs[i]);
  }
  
  for (int j = 0; j < MA->nmodels; j++){
    post_probs[j] = post_probs[j]/ norm_sum; 
    for (double  i = 0.0; i < 0.99; i += 0.01 ){
      if ( isnan(b[j].BMD_CDF.inv(i))){
        post_probs[j] = 0; // if the cdf has nan in it then it needs a 0 posterior
      }  
      
    } 
  }
  
  norm_sum = 0.0; 
  for (int i =0; i < MA->nmodels; i++){
    norm_sum += post_probs[i]; 
  }
  
  for (int i =0; i < MA->nmodels; i++){
    post_probs[i] = post_probs[i]/norm_sum; 
    res->post_probs[i] = post_probs[i];
    
    transfer_continuous_model(b[i],res->models[i]);
    transfer_mcmc_output(a[i],ma->analyses[i]); 
    res->models[i]->model = MA->models[i]; 
    res->models[i]->dist  = MA->disttype[i];
  }
  
  double range[2]; 
  /*
  bmd_range_find(res,range);
  double range_bmd = range[1] - range[0]; 
  for (int i = 0; i < res->dist_numE; i ++){
    double cbmd = double(i)/double(res->dist_numE)*range_bmd; 
    double prob = 0.0; 
    
    for (int j = 0; j < MA->nmodels; j++){
      prob += isnan(a[j].BMD_CDF.P(cbmd))?0:a[j].BMD_CDF.P(cbmd)*post_probs[j]; 
    }
    res->bmd_dist[i] = cbmd; 
    res->bmd_dist[i+res->dist_numE]  = prob;
  }
  */
  return; 
}
