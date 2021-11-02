#ifdef R_COMPILATION
//necessary things to run in R    
  #include <RcppEigen.h>
  #include <RcppGSL.h>
#else 
  #include <Eigen/Dense>
#endif

#include <vector>
#include <algorithm>
#include "polyK/polyK_setup.h"
#include "polyK/polyK.h"
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector polyk(NumericVector dose, NumericVector tumor,
                    NumericVector daysOnStudy) {
  std::vector<double> t_dose(dose.size()); 
  std::vector<int>    t_tumor(tumor.size()); 
  std::vector<int>    t_daysOnStudy(daysOnStudy.size());
  std::vector<double> temp = t_dose; 
  if (dose.size() != tumor.size() ||
      dose.size() != daysOnStudy.size()){
    stop("The variables @dose,@tumor, and @daysOnStudy need to have the same number of entries.");
  }
  // Copy over the data
  for (int i = 0; i < dose.size(); i++){
    t_dose[i] = dose[i];
    t_tumor[i] = (int) tumor[i];
    t_daysOnStudy[i] = (int)daysOnStudy[i];
    
  }
  
  PolyK::PolyKPrepareClass pkData; 
  PolyK::TDMSE_PolyK       pkTest; 
  pkData.SetupStudy(t_dose,t_tumor,t_daysOnStudy);
    Rcout << "Here 0" << std::endl; 
  pkData.prepare();
  Rcout << "Here 1" << std::endl; 
  std::sort(t_dose.begin(), t_dose.end());
  auto last = std::unique(t_dose.begin(), t_dose.end());
  t_dose.erase(last, t_dose.end());
  //
  // POLY_3_TEST
  // POLY_1PT5_TEST
  // POLY_6_TEST
    
  double t2 = pkTest.polyk_mod(pkData,
                               pkData.getNumDoseLevels(),
                               PolyK::POLY_1PT5_TEST,
                               1.0);
  Rcout << "Here 2" << std::endl; 
   Rcpp::Rcout << t2 << std::endl; 
  return dose * 2;
}



