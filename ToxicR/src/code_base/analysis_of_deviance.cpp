/* Copyright 2021  NIEHS <matt.wheeler@nih.gov>
 * 
 *
 *Permission is hereby granted, free of charge, to any person obtaining a copy of this software 
 *and associated documentation files (the "Software"), to deal in the Software without restriction, 
 *including without limitation the rights to use, copy, modify, merge, publish, distribute, 
 *sublicense, and/or sell copies of the Software, and to permit persons to whom the Software 
 *is furnished to do so, subject to the following conditions:
 *
 *The above copyright notice and this permission notice shall be included in all copies 
 *or substantial portions of the Software.
 
 *THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
 *INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
 *PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
 *HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
 *CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
 *OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 * 
 * 
 */
#include "bmd_calculate.h"
#include "bmdStruct.h"

#if defined R_COMPILATION || defined __cplusplus
  #include <cmath> 
#else
  #include <math.h>
#endif

#ifdef R_COMPILATION  
  #include <RcppEigen.h>
  #include <RcppGSL.h>
  using namespace Rcpp;
  using Rcpp::as;
#elif __cplusplus
  #include <Eigen/Dense>
#endif
  
#define NUM_PRIOR_COLS 5

#include <algorithm>
#include <vector>
#include <limits>
#include <set>


#include "lognormalTests.h"
#include "normalTests.h"
  
/* log_normal_AOD
 * 
 * 
 * 
 * 
 */
void log_normal_AOD_fits(Eigen::MatrixXd Y, Eigen::MatrixXd X, 
                         bool bSuffStat, continuous_deviance * CD){

  /////////////////////////////////////////////////////////////////////////////////////
  // Test A1
  lognormalLLTESTA1 a1Test(Y, X, bSuffStat);
  int nParms = a1Test.nParms();
  std::vector<double> fix1(nParms); for (unsigned int i = 0; i < nParms; i++) { fix1[i] = 0.0; }
  std::vector<bool> isfix1(nParms); for (unsigned int i = 0; i < nParms; i++) { isfix1[i] = false; }
  Eigen::MatrixXd a1Priors(nParms, NUM_PRIOR_COLS);
  for (unsigned int i = 0; i < nParms; i++) {
    a1Priors.row(i) << 0, 0, 1, 0, 1e8;
  }
  a1Priors.row(nParms - 1) << 0, 0, 1, -1e8, 1e8;
  
  Eigen::MatrixXd startV1(a1Test.nParms(),1); 
  for (unsigned int i = 0; i < nParms-1; i++){
    startV1(i,0) = exp(Y(i,0)); 
  }
  
  for (unsigned int i = 0; i < nParms-1; i++){
    startV1(nParms-1,0) += Y(i,2); 
  }
  startV1(nParms,0) /= double(nParms-1); 
  
  cout << startV1 << endl; 
  IDcontinuousPrior a1Init(a1Priors);
  statModel<lognormalLLTESTA1, IDcontinuousPrior> a1Model(a1Test, a1Init,
                                                          isfix1, fix1);
  optimizationResult a1Result = findMAP<lognormalLLTESTA1, IDcontinuousPrior>(&a1Model,startV1);
  cout << a1Result.max_parms << endl; 
  CD->A1 = a1Result.functionV;
  /////////////////////////////////////////////////////////////////////////////////////
  // Test A2
  lognormalLLTESTA2 a2Test(Y, X, bSuffStat);
  nParms = a2Test.nParms();
  std::vector<double> fix2(nParms); for (unsigned int i = 0; i < nParms; i++) { fix2[i] = 0.0; }
  std::vector<bool> isfix2(nParms); for (unsigned int i = 0; i < nParms; i++) { isfix2[i] = false; }
  Eigen::MatrixXd a2Priors(nParms, NUM_PRIOR_COLS);
  for (unsigned int i = 0; i < nParms; i++) {
    a2Priors.row(i) << 0, 0, 1, -1e8, 1e8;
  }
  for (unsigned int i = 0; i < nParms / 2; i++) {
    a2Priors.row(i) << 0, a1Result.max_parms(i), 1, 0, 1e8;
  }
  IDcontinuousPrior a2Init(a2Priors);

  statModel<lognormalLLTESTA2, IDcontinuousPrior> a2Model(a2Test, a2Init,
                                                          isfix2, fix2);
  optimizationResult a2Result = findMAP<lognormalLLTESTA2, IDcontinuousPrior>(&a2Model);
  CD->A2 = a2Result.functionV;
  /////////////////////////////////////////////////////////////////////////////////////
  // Test A3
  lognormalLLTESTR rTest(Y, X, bSuffStat);
  nParms = rTest.nParms();
  std::vector<double> fix4(nParms); for (unsigned int i = 0; i < nParms; i++) { fix4[i] = 0.0; }
  std::vector<bool> isfix4(nParms); for (unsigned int i = 0; i < nParms; i++) { isfix4[i] = false; }
  Eigen::MatrixXd rPriors(nParms, NUM_PRIOR_COLS);
  for (unsigned int i = 0; i < nParms; i++) {
    rPriors.row(i) << 0, 1, 1, 0, 1e8;
  }
  rPriors.row(nParms - 1) << 0, 0, 1, -1e8, 1e8;
  IDcontinuousPrior rInit(rPriors);
  statModel<lognormalLLTESTR, IDcontinuousPrior> rModel(rTest, rInit,
                                                        isfix4, fix4);
  optimizationResult rResult = findMAP<lognormalLLTESTR, IDcontinuousPrior>(&rModel);
  
  CD->A3 = rResult.functionV;
}

/* normal_AOD
 * 
 * 
 * 
 * 
 * 
 */
void normal_AOD_fits(Eigen::MatrixXd Y, Eigen::MatrixXd X, 
                     bool bSuffStat, continuous_deviance * CD){
  
  //////////////////////////////////////////////////////////////////////
  normalLLTESTA1 a1Test(Y, X, bSuffStat);
  int  nParms = a1Test.nParms();
  std::vector<double> fix1(nParms); for (unsigned int i = 0; i < nParms; i++) { fix1[i] = 0.0; }
  std::vector<bool> isfix1(nParms); for (unsigned int i = 0; i < nParms; i++) { isfix1[i] = false; }
  Eigen::MatrixXd a1Priors(nParms, NUM_PRIOR_COLS);
  for (unsigned int i = 0; i < nParms; i++) {
    a1Priors.row(i) << 0, 0, 1, -1e8, 1e8;
  }
  IDcontinuousPrior a1Init(a1Priors);
  
  statModel<normalLLTESTA1, IDcontinuousPrior> a1Model(a1Test, a1Init,
                                                       isfix1, fix1);
  optimizationResult a1Result = findMAP<normalLLTESTA1, IDcontinuousPrior>(&a1Model);
  CD->A3 = a1Result.functionV;
  ////////////////////////////////////////////////////////////////////
  normalLLTESTA2 a2Test(Y, X, bSuffStat);
  nParms = a2Test.nParms();
  std::vector<double> fix2(nParms); for (unsigned int i = 0; i < nParms; i++) { fix2[i] = 0.0; }
  std::vector<bool> isfix2(nParms); for (unsigned int i = 0; i < nParms; i++) { isfix2[i] = false; }
  Eigen::MatrixXd a2Priors(nParms, NUM_PRIOR_COLS);
  for (unsigned int i = 0; i < nParms; i++) {
    a2Priors.row(i) << 0, 0, 1, -1e8, 1e8;
  }
  for (unsigned int i = 0; i < nParms / 2; i++) {
    a2Priors.row(i) << 0, a1Result.max_parms(i), 1, -1e8, 1e8;
  }
  IDcontinuousPrior a2Init(a2Priors);
  statModel<normalLLTESTA2, IDcontinuousPrior> a2Model(a2Test, a2Init,
                                                       isfix2, fix2);
  optimizationResult a2Result = findMAP<normalLLTESTA2, IDcontinuousPrior>(&a2Model);
  CD->A2 = a2Result.functionV;
  //////////////////////////////////////////////////////////////////////
  normalLLTESTA3 a3Test(Y, X, bSuffStat);
  nParms = a3Test.nParms();
  std::vector<double> fix3(nParms); for (unsigned int i = 0; i < nParms; i++) { fix3[i] = 0.0; }
  std::vector<bool> isfix3(nParms); for (unsigned int i = 0; i < nParms; i++) { isfix3[i] = false; }
  Eigen::MatrixXd a3Priors(nParms, NUM_PRIOR_COLS);
  for (unsigned int i = 0; i < nParms; i++) {
    a3Priors.row(i) << 0, 0, 1, -1e8, 1e8;
  }
  for (unsigned int i = 0; i < nParms - 2; i++) {
    a3Priors.row(i) << 0, a2Result.max_parms(i), 1, 0, 1e8;
  }
  IDcontinuousPrior a3Init(a3Priors);
  statModel<normalLLTESTA3, IDcontinuousPrior> a3Model(a3Test, a3Init,
                                                       isfix3, fix3);
  optimizationResult a3Result = findMAP<normalLLTESTA3, IDcontinuousPrior>(&a3Model);
  CD->A3 = a3Result.functionV;
  
}