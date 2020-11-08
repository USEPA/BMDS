
// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::depends(RcppEigen)]]

#ifdef R_COMPILATION
//necessary things to run in R    
  #include <RcppEigen.h>
  #include <RcppGSL.h>
using namespace Rcpp;
#else 
  #include <Eigen/Dense>

#endif

#include <gsl/gsl_randist.h>

#include "continuous_clean_aux.h"
#include "bmdStruct.h"
#include <iostream>
#include <numeric> 


using namespace std;


// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::depends(RcppEigen)]]


/*data clean procedures*/ 
 
 std::vector<double> unique_list(Eigen::MatrixXd X){
   
   std::vector<double> uniqueX; 
   
   for (int i = 0; i < X.rows(); i++){
     bool isUnique = true; 
     for (int j = 0; j < uniqueX.size(); j++){
       if (X(i,0)==uniqueX[j]){
         isUnique = false; 
         break; 
       }
     }
     if (isUnique){
       uniqueX.push_back(X(i,0)); 
     }
   }
   
   return uniqueX; 
 }
 
 
Eigen::MatrixXd cleanSuffStat(Eigen::MatrixXd Y, Eigen::MatrixXd X, bool is_logNormal, bool use_divisor){
  double minDose = X.minCoeff(); 
  double divisor = 0; 
  int nmin = 0; 
  
  for (int i = 0; i < X.rows(); i++){
    if (X(i,0)==minDose){
      nmin++; 
      divisor += Y(i,0); 
    }
  }
  if (use_divisor){
     divisor = divisor/double(nmin); //average background dose
  }else{
     divisor = 1.0; 
  } 
  
  Y.col(0).array() = Y.col(0).array()/divisor; //divide mean
  Y.col(2).array() = Y.col(2).array()/divisor; //divide sd; 
  if (is_logNormal){
    
    Eigen::MatrixXd t1 = sqrt(log(1.0+pow(Y.col(2).array(),2)/Y.col(0).array())); 
    Eigen::MatrixXd t2 = log(Y.col(0).array())-pow(t1.array(),2)*0.5; 
    Y.col(0) = t2.array(); 
    Y.col(2) = t1.array(); 
  }
  return Y; 
  
}

double get_divisor(Eigen::MatrixXd Y, Eigen::MatrixXd X){
  // find the average of the lowest dose (e.g. background)
  double minDose = X.minCoeff(); 
  double divisor = 0; 
  int nmin = 0; 
  
  for (int i = 0; i < X.rows(); i++){
    if (X(i,0) == minDose){
      nmin++; 
      divisor += Y(i,0); 
    }
  }
  divisor = divisor/double(nmin); 
  return fabs(divisor); // return the absolute value of the divisor so we don't 
                        // flip the sign of the dose-response curve. 
}

Eigen::MatrixXd createSuffStat(Eigen::MatrixXd Y, Eigen::MatrixXd X,
                               bool is_logNormal){
  // get unique element
  std::vector<double> uniqueX = unique_list(X); 
  
  // find the average of the lowest dose (e.g. background)
  double minDose = X.minCoeff(); 
  double divisor = 0; 
  int nmin = 0; 
  for (int i = 0; i < X.rows(); i++){
    if (X(i,0)==minDose){
      nmin++; 
      divisor += Y(i,0); 
    }
  }
  divisor = divisor/double(nmin); 
  
  // build the sufficient statistics
  Eigen::MatrixXd SSTAT(uniqueX.size(),3); 
  for (int i = 0; i < uniqueX.size(); i++){
    std::vector<double> uVals; 
    for (int j = 0; j < Y.rows(); j++){
      if (X(j,0)==uniqueX[i]){
        if (is_logNormal){
          uVals.push_back(log(Y(j,0))); //geometric mean
        }else{
          uVals.push_back(Y(j,0));  // aritmetic mean
        }
      }
    }
    SSTAT(i,0) = accumulate( uVals.begin(), uVals.end(), 0.0) / uVals.size();
    SSTAT(i,1) = uVals.size();
    double sqSum = std::inner_product(uVals.begin(), uVals.end(), uVals.begin(), 0.0);
    SSTAT(i,2) =  sqSum / uVals.size() - SSTAT(i,0)*SSTAT(i,0);
    SSTAT(i,2) =  SSTAT(i,2) * (double(uVals.size())/(double(uVals.size())-1.0)); //n-1 instead of n
    SSTAT(i,2) =  pow(SSTAT(i,2),0.5); 
  }
  
  return SSTAT; 
}
 
// FIXME: CHECK IF WE ARE RESCALING THE VARIANCE PARAMETERS CORRECTLY
Eigen::MatrixXd rescale_parms(Eigen::MatrixXd parms, cont_model model,
                              double max_dose, double bkground,bool is_logNormal)
  {
    
    switch(model){
      case cont_model::hill:
        parms(0,0) *= bkground; parms(1,0) *= bkground; parms(2,0)*=max_dose; 
        if (!is_logNormal){
          if (parms.rows()==5){
            parms(4,0) += 2*log(bkground); 
          }else{
            parms(5,0) += 2*log(bkground); 
          }
        }
        break; 
      case cont_model::exp_3:
        parms(0,0) *= bkground; parms(1,0) *= 1/max_dose; 
        if (!is_logNormal){
          if (parms.rows()== 5){
            parms(4,0) += 2*log(bkground); 
          }else{
            parms(5,0) += 2*log(bkground); 
          }
        }
        break; 
      case cont_model::exp_5:
        
        parms(0,0) *= bkground; parms(1,0) *= 1/max_dose; 
        if (!is_logNormal){
            if (parms.rows()==5){
              parms(4,0) += 2*log(bkground); 
            }else{
              parms(5,0) += 2*log(bkground); 
            }
        }
        break; 
        
      case cont_model::power: 
        parms(0,0) *= bkground; parms(1,0) *= bkground; 
        parms(1,0) *= pow(1/max_dose,parms(2,0)); 
        if (!is_logNormal){
          if (parms.rows()==4){
            parms(3,0) += 2*log(bkground); 
          }else{
            parms(4,0) += 2*log(bkground); 
          }
        }
        break; 
      case cont_model::polynomial:
      
      default:
        break; 
    }
    
    return parms; 
    
  }
 
 
 Eigen::MatrixXd rescale_cov_matrix(Eigen::MatrixXd COV, 
									Eigen::MatrixXd parms, cont_model model,
									double max_dose, double bkground,
									bool is_logNormal)
  {
    Eigen::MatrixXd scaleMatrix = Eigen::MatrixXd::Identity(COV.rows(), COV.cols());
    switch(model){
      case cont_model::hill:
        scaleMatrix(0,0) = bkground; scaleMatrix(1,1) = bkground; scaleMatrix(2,2)*= max_dose; 
        COV = scaleMatrix*COV*scaleMatrix; 
        break; 
      case cont_model::exp_3:
        scaleMatrix(0,0) = bkground; scaleMatrix(1,1) = 1/max_dose;  
        COV = scaleMatrix*COV*scaleMatrix; 
        break; 
      case cont_model::exp_5:
        scaleMatrix(0,0) = bkground; scaleMatrix(1,1) = 1/max_dose; 
        COV = scaleMatrix*COV*scaleMatrix; 
        break; 
      
      case cont_model::power: 
        parms(0,0) *= bkground; parms(1,0) *= bkground*pow(1/max_dose,parms(2,0)); 
        scaleMatrix(0,0) = bkground; scaleMatrix(1,1) = bkground*pow(1/max_dose,parms(2,0));
        scaleMatrix(1,2) = bkground*parms(1,0)*log(1/max_dose)*pow(1/max_dose,parms(2,0));  
        COV = scaleMatrix*COV*scaleMatrix; 
        break; 
      case cont_model::polynomial:
        break; 
    }
    return COV; 

  }
 