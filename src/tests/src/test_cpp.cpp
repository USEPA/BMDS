// test
//
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
// #include "bmds_helper.h"
#include "priors.h"
#include "test_cpp.h"
// #include "dichotomous_entry_code.h"
// #include "continuous_entry_code.h"
// #include "bmdStruct.h"
// #include "test_cpp.h"
#include <iomanip>
#include <iostream>
#include <string>

void runOldDichoAnalysis();
void runOldContAnalysis();
void runDichoMA();
void runPythonDichoAnalysis();
void runPythonDichoMA();
void runPythonContAnalysis();
void runPythonContLoud();
void runPythonMultitumorAnalysis();
void runTestMultitumorModel();
void runPythonNestedAnalysis();
void test();
void printDichoModResult(
    struct python_dichotomous_analysis *pyAnal, struct python_dichotomous_model_result *pyRes,
    bool showResultsOverride
);
void printNestedModResult(
    struct python_nested_analysis *pyAnal, struct python_nested_result *pyRes,
    bool showResultsOverride
);
std::vector<double> getMultitumorPrior(int degree, int prior_cols);
std::vector<double> getNlogisticPrior(int ngrp, int prior_cols, bool restricted);
std::vector<double> getNctrPrior(int ngrp, int prior_cols, bool restricted);
void Nlogist_probs_test();
void Nlogist_lk_test();
void runTestMTInequalityConstraint();
void runTestMTEqualityConstraint();

bool showResultsOverride = true;

int main(void) {
  //  runOldDichoAnalysis();

  //  runOldContAnalysis();
  //  runCompleteContAnalysis();
  //  runPythonDichoAnalysis();
  //  runPythonDichoMA();
  //  runPythonContAnalysis();
  runPythonContLoud();
  //  runPythonMultitumorAnalysis();
  //  runPythonNestedAnalysis();
  //  Nlogist_probs_test();
  //  Nlogist_lk_test();
  ////  runTestMultitumorModel();
  //  runTestMTEqualityConstraint();
  return 0;
}

// struct contInputData{
//    double *D;
//    double *Y;
//    double *N;
//    double *SD;
//    bool isIncreasing;
//    bool suff_stat;
//    int numRows;
// };
//
// void getData(struct contInputData * data){
//
//   data->numRows = 9;
//   data->isIncreasing = true;
//   data->suff_stat = true;
//   data->D = malloc(data->numRows*sizeof(double));
//   data->Y = malloc(data->numRows*sizeof(double));
//   data->N = malloc(data->numRows*sizeof(double));
//   data->SD = malloc(data->numRows*sizeof(double));
//
//
//
//   double tmpD[] = {0, 0.25, 0.5, 1, 2, 4, 8, 16, 32};
//   double tmpY[] = {6.7, 6.4, 6.9, 7.5, 7.7, 9.7, 11.5, 12.9, 13.8};
//   double tmpN[] = {13, 12, 10, 13, 14, 13, 11, 13, 12};
//   double tmpSD[] = {0.360555128, 1.039230485, 0.948683298, 0.360555128, 0.748331477,
//   0.360555128, 1.658312395, 1.44222051, 3.117691454};
//
//   for (int i=0; i<data->numRows; i++){
//     data->D[i] = tmpD[i];
//     data->Y[i] = tmpY[i];
//     data->N[i] = tmpN[i];
//     data->SD[i] = tmpSD[i];
//   }
//
//   printf("output\n");
//   for (int i=0; i<data->numRows; i++){
//      printf("i:%d, %f, %f, %f, %f\n", i, data->D[i], data->Y[i], data->N[i], data->SD[i]);
//   }
//   printf("end output\n");
//
// }
//
////1=absdev, 2 = stddev, 3 = reldev, 4 = pt, 5 = extra, 6 = hybrid_extra, 7 = hybrid_added
//
// void contRuns(){
//
//
//  int numDatasets = 2;
//
//  struct contInputData* inData = malloc(numDatasets*sizeof(struct contInputData));
//
//  for (int i=0; i < numDatasets; i++){
//    getData(&inData[i]);
//  }
//
////  printf("numRows: %d\n", inData.numRows);
////  for (int i=0; i<inData.numRows; i++){
////     printf("i:%d, %f, %f, %f, %f\n", i, inData.D[i], inData.Y[i], inData.N[i], inData.SD[i]);
////  }
//
//
//    int i=0;
//    //assign everything to anal struct
//    struct continuous_analysis anal;
//    anal.isIncreasing = inData[i].isIncreasing;
//    anal.suff_stat = inData[i].suff_stat;
//    anal.doses = inData[i].D;
//    anal.Y = inData[i].Y;
//    anal.n_group = inData[i].N;
//    anal.sd = inData[i].SD;
//
//    //MODEL//////////
//    int modelType = 1;   //1 = frequentist, 2 = bayesian
//    bool restricted = false;   //only used for frequentist models
//    anal.model = power; //hill, exp_3, exp_5, power, funl, polynomial
//    anal.degree = 2;  //for polynomial only
//
//    //OPTIONS//////////
//    bool detectAdvDir = true;  //if false then need to set isIncreasing
//    anal.disttype = normal;  //normal, normal_ncv, log_normal
//    anal.alpha = 0.05;
//    anal.BMR = 1.0;
//    anal.BMD_type = 2; //1=absdev, 2 = stddev, 3 = reldev, 4 = pt, 5 = extra, 6 = hybrid_extra, 7
//    = hybrid_added
//  for (int j = 0; j<numDatasets; j++){
//
////    int i=0;
////    //assign everything to anal struct
////    struct continuous_analysis anal;
////    anal.isIncreasing = inData[i].isIncreasing;
////    anal.suff_stat = inData[i].suff_stat;
////    anal.doses = inData[i].D;
////    anal.Y = inData[i].Y;
////    anal.n_group = inData[i].N;
////    anal.sd = inData[i].SD;
////
////    //MODEL//////////
////    int modelType = 1;   //1 = frequentist, 2 = bayesian
////    bool restricted = false;   //only used for frequentist models
////    anal.model = power; //hill, exp_3, exp_5, power, funl, polynomial
////    anal.degree = 2;  //for polynomial only
////
////    //OPTIONS//////////
////    bool detectAdvDir = true;  //if false then need to set isIncreasing
////    anal.disttype = normal;  //normal, normal_ncv, log_normal
////    anal.alpha = 0.05;
////    anal.BMR = 1.0;
////    anal.BMD_type = 2; //1=absdev, 2 = stddev, 3 = reldev, 4 = pt, 5 = extra, 6 = hybrid_extra,
/// 7 = hybrid_added
//
//
//    runCompleteContAnalysis_new(&anal, inData[i].numRows, modelType, restricted, detectAdvDir);
//  }
////  for (int i = 0; i<numDatasets; i++){
////
////    //assign everything to anal struct
////    struct continuous_analysis anal;
////    anal.isIncreasing = inData[i].isIncreasing;
////    anal.suff_stat = inData[i].suff_stat;
////    anal.doses = inData[i].D;
////    anal.Y = inData[i].Y;
////    anal.n_group = inData[i].N;
////    anal.sd = inData[i].SD;
////
////    //MODEL//////////
////    int modelType = 1;   //1 = frequentist, 2 = bayesian
////    bool restricted = false;   //only used for frequentist models
////    anal.model = power; //hill, exp_3, exp_5, power, funl, polynomial
////    anal.degree = 2;  //for polynomial only
////
////    //OPTIONS//////////
////    bool detectAdvDir = true;  //if false then need to set isIncreasing
////    anal.disttype = normal;  //normal, normal_ncv, log_normal
////    anal.alpha = 0.05;
////    anal.BMR = 1.0;
////    anal.BMD_type = 2; //1=absdev, 2 = stddev, 3 = reldev, 4 = pt, 5 = extra, 6 = hybrid_extra,
/// 7 = hybrid_added
////
////
////    runCompleteContAnalysis_new(&anal, inData[i].numRows, modelType, restricted, detectAdvDir);
////  }
//
//}

///////////////////////////
// Initial values & Priors//
///////////////////////////
//
// Dichotomous Models
//
//////////BAYESIAN////////////////
//  double prBayesianLogistic[] = {1,2,0, 0, 2, 2, -20, 0, 20, 40};
//  double prBayesianDHill[] = {1,1,1,2,-1,0,-3,0.693147,2,3,3.3,0.5,-40,-40,-40,0,40,40,40,40};
//  double prBayesianGamma[] = {1,2,2,0,0.693147,0,2,0.424264,1,-18,0.2,0,18,20,10000};
//  double prBayesianLogLogistic[] = {1,1,2,0,0,0.693147,2,1,0.5,-20,-40,0,20,40,20};
//  double prBayesianLogProbit[] = {1,1,2,0,0,0.693147,2,1,0.5,-20,-40,0,20,40,20};
//  double prBayesianMulti1[] = {1,2,0,0,2,1,-20,0,20,1e6}; //degree 1
//  double prBayesianMulti2[] = {1,2,2,0,0,0,2,1,1,-20,0,0,20,1e6,1e6}; //degree 2
//  double prBayesianMulti3[] = {1,2,2,2,0,0,0,0,2,1,1,1,-20,0,0,0,20,1e6,1e6,1e6}; //degree 3
//  double prBayesianMulti4[] = {1,2,2,2,2,0,0,0,0,0,2,1,1,1,1,-20,0,0,0,0,20,1e6,1e6,1e6,1e6};
//  //degree 4 double prBayesianMulti5[] =
//  {1,2,2,2,2,2,0,0,0,0,0,0,2,1,1,1,1,1,-20,0,0,0,0,0,20,1e6,1e6,1e6,1e6,1e6}; //degree 5 double
//  prBayesianProbit[] = {1,2,0,0,2,1,-20,0,20,40}; double prBayesianQLinear[] =
//  {1,2,0,0,2,1,-20,0,20,18}; double prBayesianWeibull[] =
//  {1,2,2,0,0.424264,0,2,0.5,1.5,-20,0,0,20,40,10000};
////////////UNRESTRICTED FREQ////////////////
//  double prUFreqLogistic[] = {0,0,0,0,0,0,-18,0,18,100};
//  double prUFreqDHill[] = {0,0,0,0,0,0,0,0,0,0,0,0,-18,-18,-18,1e-8,18,18,18,18};
//  double prUFreqGamma[] = {0,0,0,0,0,0,0,0,0,-18,0.2,0,18,18,100};
//  double prUFreqLogLogistic[] = {0,0,0,0,0,0,0,0,0,-18,-18,1e-4,18,18,18};
//  double prUFreqLogProbit[] = {0,0,0,0,0,0,0,0,0,-18,-18,1e-4,18,18,18};
//  double prUFreqMulti1[] = {0,0,0,0,0,0,-18,-18,18,100}; //degree 1
//  double prUFreqMulti2[] = {0,0,0,0,0,0,0,0,0,-18,-18,-18,18,100,1e4}; //degree 2
//  double prUFreqMulti3[] = {0,0,0,0,0,0,0,0,0,0,0,0,-18,-18,-18,-18,18,100,1e4,1e4}; //degree 3
//  double prUFreqMulti4[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-18,-18,-18,-18,-18,18,100,1e4,1e4,1e4};
//  //degree 4 NOTWORKING double prUFreqMulti5[] =
//  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-18,-18,-18,-18,-18,-18,18,1e4,1e4,1e4,1e4,1e4}; //degree 5
//  NOTWORKING double prUFreqProbit[] = {0,0,0,0,0,0,-18,0,18,18}; double prUFreqQLinear[] =
//  {0,0,0,0,0,0,-18,0,18,100}; double prUFreqWeibull[] =
//  {0,0,0,0,0,0,0,0,0,-18,1e-6,1e-6,18,18,100};
////////////RESTRICTED FREQ////////////////
//  double prRFreqDHill[] = {0,0,0,0,0,0,0,0,0,0,0,0,-18,-18,-18,1,18,18,18,18};
//  double prRFreqGamma[] = {0,0,0,0,0,0,0,0,0,-18,1,0,18,18,100};
//  double prRFreqLogLogistic[] = {0,0,0,0,0,0,0,0,0,-18,-18,1,18,18,18};
//  double prRFreqLogProbit[] = {0,0,0,0,0,0,0,0,0,-18,-18,1,18,18,18};
//  double prRFreqMulti1[] = {0,0,0,0,0,0,-18,0,18,1e4}; //degree 1
//  double prRFreqMulti2[] = {0,0,0,0,0,0,0,0,0,-18,0,0,18,1e4,1e4}; //degree 2
//  double prRFreqMulti3[] = {0,0,0,0,0,0,0,0,0,0,0,0,-18,0,0,0,18,1e4,1e4,1e4}; //degree 3
//  double prRFreqMulti4[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-18,0,0,0,0,18,1e4,1e4,1e4,1e4};
//  //degree 4 double prRFreqMulti5[] =
//  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-18,0,0,0,0,0,18,1e4,1e4,1e4,1e4,1e4}; //degree 5 double
//  prRFreqWeibull[] = {0,0,0,0,0,0,0,0,0,-18,1,1e-6,18,18,100};
///////////////////////////////////////////////////////////////////////////////
////Continuous Models
///////////////////////////////////////////////////////////////////////////////
////BAYESIAN
//  double prBayesianHill[] =
//  {2,1,2,2,1,0,1,-0.69315,0.405465,0,1,2,1,0.2501,1,0,-18,0,0,-18,18,18,18,18,18}; double
//  prBayesianHillNCV[] =
//  {2,1,2,2,1,0,1,-0.69315,0.405465,0,1,2,1,0.2501,1,0,-18,0,0,-18,18,18,18,18,18}; double
//  prBayesianPower[] = {2,1,2,1,0,0,0.405465,0,1,1,0.5,1,0,-10000,0,-18,1e6,1e4,40,18}; double
//  prBayesianPowerNCV[] =
//  {2,1,2,2,1,0,0,0.405465,0,0,1,1,0.5,0.2501,1,0,-10000,0,0,-18,1e6,1e4,40,18,18};
//  //funl
//  double prBayesianPoly1[] = {2,1,1,0,0,0,1,2,1,0,-10000,-18,1e6,1e4,18};
//  double prBayesianPoly2[] = {2,1,1,1,0,0,0,0,1,2,2,1,0,-10000,-10000,-18,1e6,1e4,1e4,18}; //poly
//  2 double prBayesianPoly3[] =
//  {2,1,1,1,1,0,0,0,0,0,1,2,2,2,1,0,-10000,-10000,-10000,-18,1e6,1e4,1e4,1e4,18}; //poly 3 double
//  prBayesianPoly4[] =
//  {2,1,1,1,1,1,0,0,0,0,0,0,1,2,2,2,1,1,0,-10000,-10000,-10000,-10000,-18,1e6,1e4,1e4,1e4,1e4,18};
//  //poly 4 double prBayesianPoly5[] =
//  {2,1,1,1,1,1,1,0,0,0,0,0,0,0,1,2,2,2,2,1,1,0,-10000,-10000,-10000,-10000,-10000,-18,1e6,1e4,1e4,1e4,1e4,1e4,18};
//  //poly 5 double prBayesianPoly1NCV[] =
//  {2,1,2,1,0,0,0,0,1,2,0.2501,1,0,-10000,0,-18,1e6,1e4,18,18}; double prBayesianPoly2NCV[] =
//  {2,1,1,2,1,0,0,0,0,0,1,2,2,0.2501,1,0,-10000,-10000,0,-18,1e6,1e4,1e4,18,18}; //poly 2 double
//  prBayesianPoly3NCV[] =
//  {2,1,1,1,2,1,0,0,0,0,0,0,1,2,2,2,0.2501,1,0,-10000,-10000,-10000,0,-18,1e6,1e4,1e4,1e4,18,18};
//  //poly 3 double prBayesianPoly4NCV[] =
//  {2,1,1,1,1,2,1,0,0,0,0,0,0,0,1,2,2,2,2,0.2501,1,0,-10000,-10000,-10000,-10000,0,-18,1e6,1e4,1e4,1e4,1e4,18,18};
//  //poly 4 double prBayesianPoly5NCV[] =
//  {2,1,1,1,1,1,2,1,0,0,0,0,0,0,0,0,1,2,2,2,2,2,0.2501,1,0,-10000,-10000,-10000,-10000,-10000,0,-18,1e6,1e4,1e4,1e4,1e4,1e4,18,18};
//  //poly 5 double prBayesianExp5[] =
//  {2,2,1,2,1,0,0,0,0,0,1,1,1,0.2501,1,0,0,-20,0,-18,1e6,100,20,18,18}; double prBayesianExp5NCV[]
//  = {2,2,1,2,2,1,0,0,0,0,0,0,1,1,1,0.2501,0.5,1,0,0,-20,0,0,-18,1e6,100,20,18,18,18};
////UNRESTRICTED FREQ
//  double prUFreqHillNormal[] =
//  {0,0,0,0,0,0,0,0,1,0,1,2,1,1.2,1,-1e2,-1e2,0,1e-8,-18,1e2,1e2,5,18,18};  //normal dist double
//  prUFreqHillNormalNCV[] =
//  {0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,1,-1e8,-1e8,0,1e-8,-1e3,-1e3,1e8,1e8,30,18,1000,1000};
//  //normal dist double prUFreqHillLognormal[] =
//  {0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,1e-8,-1e8,0,1e-8,-1e3,1e8,1e8,100,100,1000};  //normal dist
//  double prUFreqPower[] = {0,0,0,0,0,0,1,0,0.1,1,0.2,1,-100,-100,0,-18,100,100,18,18};
//  double prUFreqPowerNCV[] =
//  {0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1e-8,-1e8,1e-8,-1000,-1000,1e8,1e8,100,1000,1000};
//  //funl;
//  //priors for auto detect adv dir
//  double prUFreqPoly1[] = {0,0,0,0,0,0,0,0,0,-1e6,-1e6,-18,1e6,1e6,18}; //poly 1
//  double prUFreqPoly2[] = {0,0,0,0,0,0,0,0,0,0,0,0,-1e6,-1e6,-1e6,-18,1e6,1e6,1e6,18}; //poly 2
//  double prUFreqPoly3[] =
//  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1e6,-1e6,-1e6,-1e6,-18,1e6,1e6,1e6,1e6,18}; //poly 3 double
//  prUFreqPoly4[] =
//  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1e6,-1e6,-1e6,-1e6,-1e6,-18,1e6,1e6,1e6,1e6,1e6,18};
//  //poly 4 double prUFreqPoly5[] =
//  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1e6,-1e6,-1e6,-1e6,-1e6,-1e6,-18,1e6,1e6,1e6,1e6,1e6,1e6,18};
//  //poly 5 double prUFreqPoly1NCV[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,-18,0,-18,1000,18,18,18}; //poly
//  1 double prUFreqPoly2NCV[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-18,-1e6,0,-18,1000,18,1e6,18,18};
//  //poly 2 double prUFreqPoly3NCV[] =
//  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-18,-1e6,-1e6,0,-18,1000,18,1e6,1e6,18,18}; //poly 3
//  double prUFreqPoly4NCV[] =
//  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-18,-1e6,-1e6,-1e6,0,-18,1000,18,1e6,1e6,1e6,18,18};
//  //poly 4 double prUFreqPoly5NCV[] =
//  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-18,-1e6,-1e6,-1e6,-1e6,0,-18,1000,18,1e6,1e6,1e6,1e6,18,18};
//  //poly 5
//
//
////RESTRICTED FREQ
//  double prRFreqExp5Normal[]={0,0,0,0,0,0,0,0,0,0,0.1,1,0.5,0.2,1,0,0,-20,1,-18,100,100,20,18,18};
//  double
//  prRFreqExp5NormalNCV[]={0,0,0,0,0,0,0,0,0,0,0,0,0.1,1,0.5,0.2,0.5,1,0,0,-20,1,0,-18,100,100,20,18,18,18};
//  double
//  prRFreqExp5Lognormal[]={0,0,0,0,0,0,0,0,1,0,0.1,1,1,0.2,1,-1000,0,-20,1,-18,1000,100,20,18,18};
//  double prRFreqHillNormal[] =
//  {0,0,0,0,0,0,0,0,0,0,1,2,1,1.2,1,-100,-100,0,1,-18,100,100,5,18,18};  //normal dist double
//  prRFreqHillNormalNCV[] =
//  {0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,1,-1e8,-1000,0,1,-1e3,-1e3,1e8,1000,30,18,1000,1000};
//  //normal dist double prRFreqHillLognormal[] =
//  {0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,1e-8,-1e8,0,1,-1e3,1e8,1e8,100,100,1000};  //normal dist double
//  prRFreqPower[] = {0,0,0,0,0,0,0,0,0.1,1,0.2,1,-100,-100,1,-18,100,100,18,18};  //SEG FAULT
//  double prRFreqPowerNCV[] = {0,0,0,0,0,0,1,0,1,1,0.2,1,-100,-10000,0,-18,100,10000,18,18};  //SEG
//  FAULT
//  //funl;
//  //priors for auto detect adv dir
//  double prRFreqPoly1[] = {0,0,0,0,0,0,5,5,1,-1000,-18,-18,1000,18,18}; //poly 1
//  double prRFreqPoly2[] = {0,0,0,0,0,0,0,0,5,5,5,1,-1e6,-1e6,-1e6,-18,1e6,1e6,1e6,18}; //poly 2
//  double prRFreqPoly3[] =
//  {0,0,0,0,0,0,0,0,0,0,5,5,5,5,1,-1e6,-1e6,-1e6,-1e6,-18,1e6,1e6,1e6,1e6,18}; //poly 3 double
//  prRFreqPoly4[] =
//  {0,0,0,0,0,0,0,0,0,0,0,0,5,5,5,5,5,1,-1e6,-1e6,-1e6,-1e6,-1e6,-18,1e6,1e6,1e6,1e6,1e6,18};
//  //poly 4 double prRFreqPoly5[] =
//  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,5,5,5,5,5,5,1,-1e6,-1e6,-1e6,-1e6,-1e6,-1e6,-18,1e6,1e6,1e6,1e6,1e6,1e6,18};
//  //poly 5 double prRFreqPoly1NCV[] =
//  {0,0,0,0,0,0,0,0,1,1,1,1,-1e8,-1e8,1000,-1e8,1e8,1e8,1000,1e8}; //poly 1 double
//  prRFreqPoly2NCV[] = {0,0,0,0,0,0,0,0,0,0,5,5,5,1,1,0,-18,-1e6,0,-18,1000,18,1e6,18,18}; //poly 2
//  double prRFreqPoly3NCV[] =
//  {0,0,0,0,0,0,0,0,0,0,0,0,5,5,5,5,1,1,-1000,-10000,-10000,-10000,0,-18,1000,10000,10000,10000,100,18};
//  //poly 3 double prRFreqPoly4NCV[] =
//  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,-1e8,-1e8,-1e8,-1e8,-1e8,-1000,-1e8,1e8,1e8,1e8,1e8,1e8,1000,1e8};
//  //poly 4 double prRFreqPoly5NCV[] =
//  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,-1e8,-1e8,-1e8,-1e8,-1e8,-1e8,-1000,-1e8,1e8,1e8,1e8,1e8,1e8,1e8,1000,1e8};
//  //poly 5
//  //priors for adv dir down
//  double prRFreqPoly1Down[] = {0,0,0,0,0,0,1,1,1,-1e8,-1e8,-1e8,1e8,0,0}; //poly 1
//  double prRFreqPoly2Down[] = {0,0,0,0,0,0,0,0,1,1,1,1,-1e8,-1e8,-1e8,-1e8,1e8,0,0,0}; //poly 2
//  double prRFreqPoly3Down[] =
//  {0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,-1e8,-1e8,-1e8,-1e8,-1e8,1e8,0,0,0,0}; //poly 3 double
//  prRFreqPoly4Down[] =
//  {0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,-1e8,-1e8,-1e8,-1e8,-1e8,-1e8,1e8,0,0,0,0,0}; //poly 4
//  double prRFreqPoly5Down[] =
//  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,-1e8,-1e8,-1e8,-1e8,-1e8,-1e8,-1e8,1e8,0,0,0,0,0,0};
//  //poly 5 double prRFreqPoly1NCVDown[] =
//  {0,0,0,0,0,0,0,0,1,1,1,1,-1e8,-1e8,1000,-1e8,1e8,0,1000,0}; //poly 1 double
//  prRFreqPoly2NCVDown[] =
//  {0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,-1e8,-1e8,-1e8,-1000,-1e8,1e8,0,0,1000,0}; //poly 2 double
//  prRFreqPoly3NCVDown[] =
//  {0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,-1e8,-1e8,-1e8,-1e8,-1000,-1e8,1e8,0,0,0,1000,0}; //poly 3
//  double prRFreqPoly4NCVDown[] =
//  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,-1e8,-1e8,-1e8,-1e8,-1e8,-1000,-1e8,1e8,0,0,0,0,1000,0};
//  //poly 4 double prRFreqPoly5NCVDown[] =
//  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,-1e8,-1e8,-1e8,-1e8,-1e8,-1e8,-1000,-1e8,1e8,0,0,0,0,0,1000,0};
//  //poly 5
//
//
//
//  //priors for adv dir up
//  double prRFreqPoly1Up[] = {0,0,0,0,0,0,1,1,1,-1e8,0,0,1e8,1e8,1e8}; //poly 1
//  double prRFreqPoly2Up[] = {0,0,0,0,0,0,0,0,1,1,1,1,-1e8,0,0,0,1e8,1e8,1e8,1e8}; //poly 2
//  double prRFreqPoly3Up[] = {0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,-1e8,0,0,0,0,1e8,1e8,1e8,1e8,1e8};
//  //poly 3 double prRFreqPoly4Up[] =
//  {0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,-1e8,0,0,0,0,0,1e8,1e8,1e8,1e8,1e8,1e8}; //poly 4 double
//  prRFreqPoly5Up[] =
//  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,-1e8,0,0,0,0,0,0,1e8,1e8,1e8,1e8,1e8,1e8,1e8}; //poly
//  5
//
//
////old
//  double prRFreqPoly1NCVUp[] = {0,0,0,0,0,0,0,0,1,1,1,1,-1e8,0,-1000,0,1e8,1e8,1000,1e8}; //poly 1
////new
//  double prRFreqPoly2NCVUp[] =
//  {0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,-1e8,0,0,-1000,-1000,1e8,1e8,1e8,1000,1000}; //poly 2 double
//  prRFreqPoly3NCVUp[] =
//  {0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,-1e8,0,0,0,-1000,-1000,1e8,1e8,1e8,1e8,1000,1000}; //poly 3
//  double prRFreqPoly4NCVUp[] =
//  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,-1e8,0,0,0,0,-1000,-1000,1e8,1e8,1e8,1e8,1e8,1000,1000};
//  //poly 4 double prRFreqPoly5NCVUp[] =
//  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,-1e8,0,0,0,0,0,-1000,-1000,1e8,1e8,1e8,1e8,1e8,1e8,1000,1000};
//  //poly 5

void runOldDichoAnalysis() {
  printf("Running dichotomous analysis\n");

  ///////////////////////////////
  // USER INPUT
  ///////////////////////////////

  enum dich_model model = d_weibull;  // d_hill =1, d_gamma=2,d_logistic=3, d_loglogistic=4,
                                      // d_logprobit=5, d_multistage=6,d_probit=7,
                                      // d_qlinear=8,d_weibull=9
  int modelType = 1;                  // 1 = frequentist, 2 = bayesian
  bool restricted = true;             // only used for frequentist models
  int BMD_type = 1;                   // 1 = extra ; added otherwise
  int degree = 3;                     // for multistage only
  double BMR = 0.1;
  double alpha = 0.05;
  bool countAllParmsOnBoundary = true;
  ///////////////////////////////
  // dicho data - dose, N, incidence
  ///////////////////////////////
  //  //Dichotomous.dax Effect 1
  //  double D[] = {0,50, 100, 150, 200};
  //  double Y[] = {0, 5, 30, 65, 90};
  //  double N[] = {100, 100, 100, 100, 100};

  //  //Dichotomous.dax Effect 2
  //  double D[] = {0,50, 100, 150, 200};
  //  double Y[] = {5, 10, 33, 67, 93};
  //  double N[] = {100, 100, 100, 100, 100};

  //  double D[] = {0,10, 20, 40};
  //  double Y[] = {0, 0, 2, 4};
  //  double N[] = {10, 10, 10, 10};

  //  double D[] = {12, 15, 18, 21};
  //  double Y[] = {1, 2, 3, 4};
  //  double N[] = {5, 6, 7, 8};

  //  double D[] = {0, 10, 50, 150, 400};
  //  double Y[] = {0, 0, 1, 4, 11};
  //  double N[] = {20, 20, 20, 20, 20};

  // Allen test data
  //  double D[] = {0, 10, 50, 150, 400};
  //  double Y[] = {0, 0, 1, 4, 11};
  //  double N[] = {20, 20, 20, 20, 20};

  // test data for extra/added risk
  // double D[] = {5,50, 100, 150, 200};
  // double Y[] = {2, 5, 30, 65, 90};
  // double N[] = {100, 100, 100, 100, 100};

  //  double D[] = {0, 0.078, 0.195};
  //  double Y[] = {0, 0, 28};
  //  double N[] = {126, 25, 119};

  //  double D[] = {0,2364.7,4973.5};
  //  double Y[] = {22,24,33};
  //  double N[] = {50,47,47};

  //  double D[] = {0, 475.1, 992.4};
  //  double Y[] = {5, 27, 40};
  //  double N[] = {50, 47, 47};

  //  //GSL Gamma crash
  //  double D[] = {0,10, 50, 150, 800};
  //  double Y[] = {0, 0, 1, 4, 11};
  //  double N[] = {20, 20, 20, 20, 20};

  // D1
  //  double D[] = {0, 0.078, 0.195};
  //  double Y[] = {0, 0, 28};
  //  double N[] = {126, 25, 119};

  // D2
  //  double D[] = {0, 18.4, 27.8};
  //  double Y[] = {2, 43, 40};
  //  double N[] = {50, 49, 45};

  // D60
  //  double D[] = {0, 0.011, 0.057, 1.3, 5.6};
  //  double Y[] = {30, 26, 17, 27, 42};
  //  double N[] = {127, 63, 64, 64, 64};

  // D80
  //  double D[] = {0, 4.79, 9.57};
  //  double Y[] = {0, 7, 8};
  //  double N[] = {10, 15, 24};

  //  //D103
  //  double D[] = {0, 786.8, 846, 925.7};
  //  double Y[] = {3, 8, 14, 23};
  //  double N[] = {50, 50, 50, 50};

  //  user submitted
  //  double D[] = {0, 11, 30, 100, 356};
  //  double Y[] = {2, 10, 13, 15, 15};
  //  double N[] = {14, 15, 15, 15, 15};

  // Bruce F Rat 2-yr C-cell Adenoma
  //  double D[] = {0, 18.2, 39.3, 74.3};
  //  double Y[] = {5, 13, 13, 8};
  //  double N[] = {50, 50, 49, 50};

  // Nasal Lesions - 2-EHA
  double D[] = {0, 10, 30, 100};
  double Y[] = {0, 0, 8, 20};
  double N[] = {20, 20, 20, 20};

  /////////////////////////////////////////////////
  ////END USER INPUT
  ////////////////////////////////////////////////////

  struct dichotomous_analysis anal;

  int numDataRows = sizeof(D) / sizeof(D[0]);

  // check data array sizes for consistency
  size_t numElementsY = sizeof(Y) / sizeof(Y[0]);
  size_t numElementsN = sizeof(N) / sizeof(N[0]);
  if (numDataRows != numElementsY || numElementsY != numElementsN) {
    printf("Number of data elements are not consistent\nExiting Code\n");
    exit(-1);
  }

  // priors defined columnwise
  int prCols = 5;

  // define priors/parameter constraints
  int numParms;
  printf("model = %d\n", model);
  switch (model) {
    case d_hill:
      numParms = 4;
      break;
    case d_gamma:
      numParms = 3;
      break;
    case d_logistic:
      numParms = 2;
      break;
    case d_loglogistic:
      numParms = 3;
      break;
    case d_logprobit:
      numParms = 3;
      break;
    case d_multistage:
      // numParms = 2 + degree;
      numParms = 1 + degree;
      break;
    case d_probit:
      numParms = 2;
      break;
    case d_qlinear:
      numParms = 2;
      break;
    case d_weibull:
      numParms = 3;
      break;
    default:
      printf("error in numParms\n");
      return;
  }

  if (modelType == 1) {
    // frequentist
    if (restricted) {
      switch (model) {
        case d_hill:
          anal.model = d_hill;
          anal.prior = prRFreqDHill;
          break;
        case d_gamma:
          anal.model = d_gamma;
          anal.prior = prRFreqGamma;
          break;
        case d_logistic:
          printf("error with restricted logistic model\n");
          return;
          break;
        case d_loglogistic:
          anal.model = d_loglogistic;
          anal.prior = prRFreqLogLogistic;
          break;
        case d_logprobit:
          anal.model = d_logprobit;
          anal.prior = prRFreqLogProbit;
          break;
        case d_multistage:
          anal.model = d_multistage;
          if (degree == 1) {
            anal.prior = prRFreqMulti1;
          } else if (degree == 2) {
            anal.prior = prRFreqMulti2;
          } else if (degree == 3) {
            anal.prior = prRFreqMulti3;
          } else if (degree == 4) {
            anal.prior = prRFreqMulti4;
          } else if (degree == 5) {
            anal.prior = prRFreqMulti5;
          }
          break;
        case d_probit:
          printf("error with restricted probit model\n");
          return;
          break;
        case d_qlinear:
          printf("error with restricted QLinear model\n");
          return;
          break;
        case d_weibull:
          anal.model = d_weibull;
          anal.prior = prRFreqWeibull;
          break;
        default:
          printf("error with restricted models\n");
          return;
      }
    } else {
      // unrestricted
      switch (model) {
        case d_hill:
          anal.model = d_hill;
          anal.prior = prUFreqDHill;
          break;
        case d_gamma:
          anal.model = d_gamma;
          anal.prior = prUFreqGamma;
          break;
        case d_logistic:
          anal.model = d_logistic;
          anal.prior = prUFreqLogistic;
          break;
        case d_loglogistic:
          anal.model = d_loglogistic;
          anal.prior = prUFreqLogLogistic;
          break;
        case d_logprobit:
          anal.model = d_logprobit;
          anal.prior = prUFreqLogProbit;
          break;
        case d_multistage:
          anal.model = d_multistage;
          if (degree == 1) {
            anal.prior = prUFreqMulti1;
          } else if (degree == 2) {
            anal.prior = prUFreqMulti2;
          } else if (degree == 3) {
            anal.prior = prUFreqMulti3;
          } else if (degree == 4) {
            anal.prior = prUFreqMulti4;
          } else if (degree == 5) {
            anal.prior = prUFreqMulti5;
          }
          break;
        case d_probit:
          anal.model = d_probit;
          anal.prior = prUFreqProbit;
          break;
        case d_qlinear:
          anal.model = d_qlinear;
          anal.prior = prUFreqQLinear;
          break;
        case d_weibull:
          anal.model = d_weibull;
          anal.prior = prUFreqWeibull;
          break;
        default:
          printf("error with restricted models\n");
          return;
      }
    }
  } else {
    // bayesian
    switch (model) {
      case d_hill:
        anal.model = d_hill;
        anal.prior = prBayesianDHill;
        break;
      case d_gamma:
        anal.model = d_gamma;
        anal.prior = prBayesianGamma;
        break;
      case d_logistic:
        anal.model = d_logistic;
        anal.prior = prBayesianLogistic;
        return;
        break;
      case d_loglogistic:
        anal.model = d_loglogistic;
        anal.prior = prBayesianLogLogistic;
        break;
      case d_logprobit:
        anal.model = d_logprobit;
        anal.prior = prBayesianLogProbit;
        break;
      case d_multistage:
        anal.model = d_multistage;
        if (degree == 1) {
          anal.prior = prBayesianMulti1;
        } else if (degree == 2) {
          anal.prior = prBayesianMulti2;
        } else if (degree == 3) {
          anal.prior = prBayesianMulti3;
        } else if (degree == 4) {
          anal.prior = prBayesianMulti4;
        } else if (degree == 5) {
          anal.prior = prBayesianMulti5;
        }
        break;
      case d_probit:
        anal.model = d_probit;
        anal.prior = prBayesianProbit;
        // printf("error with restricted probit model\n");
        // return;
        break;
      case d_qlinear:
        anal.model = d_qlinear;
        anal.prior = prBayesianQLinear;
        // printf("error with restricted QLinear model\n");
        // return;
        break;
      case d_weibull:
        anal.model = d_weibull;
        anal.prior = prBayesianWeibull;
        break;
      default:
        printf("error with restricted models\n");
        return;
    }
  }

  //////////BAYESIAN////////////////
  //  anal.model = d_logistic;
  //  double pr[] = {1,2,0, 0.1, 2, 1, -20, 1e-12, 20, 100};
  //  anal.model = d_hill;
  //  double pr[] = {1,1,1,2,-1,0,-3,0.693147,2,3,3.3,0.5,-40,-40,-40,1e-8,40,40,40,40};
  //  anal.model = d_gamma;
  //  double pr[] = {1,2,2,0,0.693147,0,2,0.424264,1,-18,0.2,1e-4,18,20,100};
  //  anal.model = d_loglogistic;
  //  double pr[] = {1,1,2,0,0,0.693147,2,1,0.5,-20,-40,1e-4,20,40,20};
  //  anal.model = d_logprobit;
  //  double pr[] = {1,1,2,0,0,0.693147,2,1,0.5,-20,-8,1e-4,20,8,40};
  //  anal.model = d_multistage;
  //  double pr[] = {1,2,0,0,2,0.5,-20,1e-4,20,100}; //degree 1
  //  double pr[] = {1,2,2,0,0,0,2,0.5,1,-20,1e-4,1e-4,20,100,1e6}; //degree 2
  //  double pr[] = {1,2,2,2,0,0,0,0,2,0.5,1,1,-20,1e-4,1e-4,1e-4,20,100,1e6,1e6}; //degree 3
  //  double pr[] = {1,2,2,2,2,0,0,0,0,0,2,0.5,1,1,1,-20,1e-4,1e-4,1e-4,1e-4,20,100,1e6,1e6,1e6};
  //  //degree 4 anal.model = d_probit; double pr[] = {1,2,0,0.1,2,1,-8,0,8,40}; anal.model =
  //  d_qlinear; double pr[] = {1,2,0,0.5,2,1,-20,0,20,100}; anal.model = d_weibull; double pr[] =
  //  {1,2,2,0,0.693147,0,2,0.424264,1,-20,1e-4,1e-4,20,18,20};

  //////////UNRESTRICTED FREQ////////////////
  //  anal.model = d_logistic;
  //  double pr[] = {0,0,0,0,0,0,-18,0,18,100};
  //  anal.model = d_hill;
  //  double pr[] = {0,0,0,0,0,0,0,0,0,0,0,0,-18,-18,-18,1e-8,18,18,18,18};
  //  anal.model = d_gamma;
  //  double pr[] = {0,0,0,0,0,0,0,0,0,-18,0.2,0,18,18,100};
  //  anal.model = d_loglogistic;
  //  double pr[] = {0,0,0,0,0,0,0,0,0,-18,-18,1e-4,18,18,18};
  //  anal.model = d_logprobit;
  //  double pr[] = {0,0,0,0,0,0,0,0,0,-18,-18,1e-4,18,18,18};
  //  anal.model = d_multistage;
  //  double pr[] = {0,0,0,0,0,0,0,0,0,-18,-18,-18,18,100,1e4};
  //  double pr[] = {0,0,0,0,0,0,-18,-18,18,100}; //degree 1
  //  double pr[] = {0,0,0,0,0,0,0,0,0,-18,-18,-18,18,100,1e4}; //degree 2
  //  double pr[] = {0,0,0,0,0,0,0,0,0,0,0,0,-18,-18,-18,-18,18,100,1e4,1e4}; //degree 3
  //  double pr[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-18,-18,-18,-18,-18,18,100,1e4,1e4,1e4}; //degree
  //  4 NOTWORKING anal.model = d_probit; double pr[] = {0,0,0,0,0,0,-18,0,18,18}; anal.model =
  //  d_qlinear; double pr[] = {0,0,0,0,0,0,-18,0,18,100}; anal.model = d_weibull; double pr[] =
  //  {0,0,0,0,0,0,0,0,0,-18,1e-6,1e-6,18,18,100};

  //////////RESTRICTED FREQ////////////////
  //  anal.model = d_hill;
  //  double pr[] = {0,0,0,0,0,0,0,0,0,0,0,0,-18,-18,-18,1,18,18,18,18};
  //  anal.model = d_gamma;
  //  double pr[] = {0,0,0,0,0,0,0,0,0,-18,1,0,18,18,100};
  //  anal.model = d_loglogistic;
  //  double pr[] = {0,0,0,0,0,0,0,0,0,-18,-18,1,18,18,18};
  //  anal.model = d_logprobit;
  //  double pr[] = {0,0,0,0,0,0,0,0,0,-18,-18,1,18,18,18};
  //  anal.model = d_multistage;
  //  double pr[] = {0,0,0,0,0,0,-18,0,18,100}; //degree 1
  //  double pr[] = {0,0,0,0,0,0,0,0,0,-18,0,0,18,100,1e4}; //degree 2
  //  double pr[] = {0,0,0,0,0,0,0,0,0,0,0,0,-18,0,0,0,18,100,1e4,1e4}; //degree 3
  // NOTWORKING  double pr[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-18,0,0,0,0,18,100,1e4,1e4,1e4};
  // //degree 4
  //  anal.model = d_weibull;
  //  double pr[] = {0,0,0,0,0,0,0,0,0,-18,1,1e-6,18,18,100};

  // parms array declared
  // int numParms = sizeof(pr)/sizeof(pr[0])/prCols;
  // double parms[numParms];
  double *parms = new double[numParms];

  // declare analysis
  anal.BMD_type = BMD_type;
  anal.BMR = BMR;
  anal.alpha = alpha;
  anal.parms = numParms;
  anal.Y = Y;
  anal.n_group = N;
  anal.doses = D;
  // anal.prior = pr;
  anal.prior_cols = prCols;
  anal.n = numDataRows;
  anal.degree = degree;

  //  printf("numParms= %d\n",numParms);
  //  printf("prCols= %d\n",prCols);

  //  if (anal.model == d_multistage){
  //    anal.degree = anal.parms - 1;
  //    printf("running multistage degree: %d\n", anal.degree);
  //  }

  // struct dichotomous_model_result *res;
  // res = new_dichotomous_model_result(anal.model, anal.parms, 200);
  struct dichotomous_model_result res;
  res.model = anal.model;
  res.parms = parms;
  res.dist_numE = 200;
  res.nparms = anal.parms;

  // double cov[numParms*numParms];
  // double bmd_dist[res.dist_numE*2];
  double *cov = new double[numParms * numParms];
  double *bmd_dist = new double[res.dist_numE * 2];

  res.cov = cov;
  res.bmd_dist = bmd_dist;

  // struct dichotomous_PGOF_result gofRes;
  struct dichotomous_GOF gof;

  struct BMDS_results bmdsRes;

  // bool* bounded = new bool[anal.parms];
  //  double* stdErr = new double[anal.parms];
  // double* lowerConf = new double[anal.parms];
  // double* upperConf = new double[anal.parms];

  // set all parms as unbounded initially
  for (int i = 0; i < anal.parms; i++) {
    // bounded[i] = false;
    bmdsRes.bounded.push_back(false);
    // stdErr[i] = -9999.0;
    bmdsRes.stdErr.push_back(BMDS_MISSING);
    bmdsRes.lowerConf.push_back(BMDS_MISSING);
    bmdsRes.upperConf.push_back(BMDS_MISSING);
    // lowerConf[i] = -9999.0;
    // upperConf[i] = -9999.0;
  }
  // bmdsRes.bounded = bounded;
  //  bmdsRes.stdErr = stdErr;
  // bmdsRes.lowerConf = lowerConf;
  // bmdsRes.upperConf = upperConf;
  bmdsRes.BMD = -9999.0;
  bmdsRes.BMDU = -9999.0;
  bmdsRes.BMDL = -9999.0;
  bmdsRes.AIC = -9999.0;

  struct dicho_AOD aod;
  double A1;
  int N1;
  double A2;
  int N2;
  double fittedLL;
  int NFit;
  double devFit;
  double devRed;
  int dfFit;
  int dfRed;
  int pvFit;
  int pvRed;
  aod.fullLL = A1;
  aod.nFull = N1;
  aod.redLL = A2;
  aod.nRed = N2;
  aod.fittedLL = fittedLL;
  aod.nFit = NFit;
  aod.devFit = devFit;
  aod.devRed = devRed;
  aod.dfFit = dfFit;
  aod.dfRed = dfRed;
  aod.pvFit = pvFit;
  aod.pvRed = pvRed;

  runBMDSDichoAnalysis(&anal, &res, &gof, &bmdsRes, &aod, &countAllParmsOnBoundary);

  printf("tlink bmdsRes.validResult = %s\n", bmdsRes.validResult ? "valid" : "invalid");
  if (bmdsRes.validResult || showResultsOverride) {
    printf("\nBenchmark Dose\n");
    printf("max: %f\n", res.max);
    printf("BMD: %f\n", bmdsRes.BMD);
    printf("BMDL: %f\n", bmdsRes.BMDL);
    printf("BMDU: %f\n", bmdsRes.BMDU);
    printf("AIC: %f\n", bmdsRes.AIC);
    printf("LPP: %f\n", bmdsRes.BIC_equiv);
    printf("P-value: %f\n", gof.p_value);
    printf("DOF: %f\n", gof.df);
    printf("Chi^2: %f\n", bmdsRes.chisq);

    //  calcDichoAIC(&anal, &res, &bmdsRes);
    printf("\nModel Parameters\n");
    printf("# of parms: %d\n", anal.parms);
    printf("parm, estimate, bounded, std.err., lower conf, upper conf\n");
    for (int i = 0; i < anal.parms; i++) {
      printf(
          "%d, %.10f, %s, %f, %f, %f\n", i, res.parms[i], bmdsRes.bounded[i] ? "true" : "false",
          bmdsRes.stdErr[i], bmdsRes.lowerConf[i], bmdsRes.upperConf[i]
      );
    }

    printf("\ncov matrix\n");
    for (int i = 0; i < anal.parms * anal.parms; i++) {
      printf("%d, %f\n", i, cov[i]);
    }

    printf("\nGoodness of Fit\n");
    printf("Dose, EstProb, Expected, Observed, Size, ScaledRes\n");
    for (int i = 0; i < gof.n; i++) {
      printf(
          "%f, %f, %f, %f, %f, %f\n", anal.doses[i], gof.expected[i] / anal.n_group[i],
          gof.expected[i], anal.Y[i], anal.n_group[i], gof.residual[i]
      );
    }
    printf("\nError bars\n");
    for (int i = 0; i < gof.n; i++) {
      printf("%f, %f\n", gof.ebLower[i], gof.ebUpper[i]);
    }

    printf("\nAnalysis of Deviance\n");
    printf("  Model,   LL,    #parms,   deviance,   test DF,  pval\n");
    printf("Full Model,  %f,  %d,  -,  -,  NA\n", aod.fullLL, aod.nFull);
    printf(
        "Fitted Model,  %f,  %d,  %f,  %d,  %f\n", aod.fittedLL, aod.nFit, aod.devFit, aod.dfFit,
        aod.pvFit
    );
    printf(
        "Reduced Model,  %f,  %d,  %f,  %d,  %f\n", aod.redLL, aod.nRed, aod.devRed, aod.dfRed,
        aod.pvRed
    );

    printf("\nBMD Dist:\n");
    for (int i = 0; i < res.dist_numE; i++) {
      printf("i:%d, perc:%f, dist:%f\n", i, res.bmd_dist[i + res.dist_numE], res.bmd_dist[i]);
    }
  } else {
    printf("\nModel was not run\n");
  }
}

void runPythonDichoAnalysis() {
  printf("Running dichotomous analysis\n");

  ///////////////////////////////
  // USER INPUT
  ///////////////////////////////

  enum dich_model model = d_multistage;  // d_hill =1, d_gamma=2,d_logistic=3, d_loglogistic=4,
                                         // d_logprobit=5, d_multistage=6,d_probit=7,
                                         // d_qlinear=8,d_weibull=9
  int modelType = 1;                     // 1 = frequentist, 2 = bayesian
  bool restricted = true;                // only used for frequentist models
  int BMD_type = 1;                      // 1 = extra ; added otherwise
  int degree = 2;                        // for multistage only
  double BMR = 0.1;
  double alpha = 0.05;
  double countAllParmsOnBoundary = false;
  ///////////////////////////////
  // dicho data - dose, N, incidence
  ///////////////////////////////
  // Dichotomous.dax Effect 1
  double D[] = {0, 50, 100, 150, 200};
  double Y[] = {0, 5, 30, 65, 90};
  double N[] = {100, 100, 100, 100, 100};

  //  //Dichotomous.dax Effect 2
  //  double D[] = {0,50, 100, 150, 200};
  //  double Y[] = {5, 10, 33, 67, 93};
  //  double N[] = {100, 100, 100, 100, 100};

  //  double D[] = {0,10, 20, 40};
  //  double Y[] = {0, 0, 2, 4};
  //  double N[] = {10, 10, 10, 10};

  //  double D[] = {12, 15, 18, 21};
  //  double Y[] = {1, 2, 3, 4};
  //  double N[] = {5, 6, 7, 8};

  //  double D[] = {0, 10, 50, 150, 400};
  //  double Y[] = {0, 0, 1, 4, 11};
  //  double N[] = {20, 20, 20, 20, 20};

  // Allen test data
  //  double D[] = {0, 10, 50, 150, 400};
  //  double Y[] = {0, 0, 1, 4, 11};
  //  double N[] = {20, 20, 20, 20, 20};

  // test data for extra/added risk
  // double D[] = {5,50, 100, 150, 200};
  // double Y[] = {2, 5, 30, 65, 90};
  // double N[] = {100, 100, 100, 100, 100};

  //  double D[] = {0, 0.078, 0.195};
  //  double Y[] = {0, 0, 28};
  //  double N[] = {126, 25, 119};

  //  double D[] = {0,2364.7,4973.5};
  //  double Y[] = {22,24,33};
  //  double N[] = {50,47,47};

  //  double D[] = {0, 475.1, 992.4};
  //  double Y[] = {5, 27, 40};
  //  double N[] = {50, 47, 47};

  //  //GSL Gamma crash
  //  double D[] = {0,10, 50, 150, 800};
  //  double Y[] = {0, 0, 1, 4, 11};
  //  double N[] = {20, 20, 20, 20, 20};

  // D1
  //  double D[] = {0, 0.078, 0.195};
  //  double Y[] = {0, 0, 28};
  //  double N[] = {126, 25, 119};

  // D2
  //  double D[] = {0, 18.4, 27.8};
  //  double Y[] = {2, 43, 40};
  //  double N[] = {50, 49, 45};

  // D60
  //  double D[] = {0, 0.011, 0.057, 1.3, 5.6};
  //  double Y[] = {30, 26, 17, 27, 42};
  //  double N[] = {127, 63, 64, 64, 64};

  // D80
  //  double D[] = {0, 4.79, 9.57};
  //  double Y[] = {0, 7, 8};
  //  double N[] = {10, 15, 24};

  //  //D103
  //  double D[] = {0, 786.8, 846, 925.7};
  //  double Y[] = {3, 8, 14, 23};
  //  double N[] = {50, 50, 50, 50};

  //  user submitted
  //  double D[] = {0, 11, 30, 100, 356};
  //  double Y[] = {2, 10, 13, 15, 15};
  //  double N[] = {14, 15, 15, 15, 15};

  // Bruce F Rat 2-yr C-cell Adenoma
  //  double D[] = {0, 18.2, 39.3, 74.3};
  //  double Y[] = {5, 13, 13, 8};
  //  double N[] = {50, 50, 49, 50};

  // Nasal Lesions - 2-EHA
  //  double D[] = {0, 10, 30, 100};
  //  double Y[] = {0, 0, 8, 20};
  //  double N[] = {20, 20, 20, 20};

  /////////////////////////////////////////////////
  ////END USER INPUT
  ////////////////////////////////////////////////////

  // struct dichotomous_analysis anal;
  struct python_dichotomous_analysis anal;
  anal.countAllParmsOnBoundary = countAllParmsOnBoundary;

  int numDataRows = sizeof(D) / sizeof(D[0]);

  // check data array sizes for consistency
  size_t numElementsY = sizeof(Y) / sizeof(Y[0]);
  size_t numElementsN = sizeof(N) / sizeof(N[0]);
  if (numDataRows != numElementsY || numElementsY != numElementsN) {
    printf("Number of data elements are not consistent\nExiting Code\n");
    exit(-1);
  }

  // priors defined columnwise
  int prCols = 5;

  // define priors/parameter constraints
  int numParms;
  printf("model = %d\n", model);
  switch (model) {
    case d_hill:
      numParms = 4;
      break;
    case d_gamma:
      numParms = 3;
      break;
    case d_logistic:
      numParms = 2;
      break;
    case d_loglogistic:
      numParms = 3;
      break;
    case d_logprobit:
      numParms = 3;
      break;
    case d_multistage:
      // numParms = 2 + degree;
      numParms = 1 + degree;
      break;
    case d_probit:
      numParms = 2;
      break;
    case d_qlinear:
      numParms = 2;
      break;
    case d_weibull:
      numParms = 3;
      break;
    default:
      printf("error in numParms\n");
      return;
  }

  double *prior;
  if (modelType == 1) {
    // frequentist
    if (restricted) {
      switch (model) {
        case d_hill:
          anal.model = d_hill;
          prior = prRFreqDHill;
          break;
        case d_gamma:
          anal.model = d_gamma;
          prior = prRFreqGamma;
          break;
        case d_logistic:
          printf("error with restricted logistic model\n");
          return;
          break;
        case d_loglogistic:
          anal.model = d_loglogistic;
          prior = prRFreqLogLogistic;
          break;
        case d_logprobit:
          anal.model = d_logprobit;
          prior = prRFreqLogProbit;
          break;
        case d_multistage:
          anal.model = d_multistage;
          if (degree == 1) {
            prior = prRFreqMulti1;
          } else if (degree == 2) {
            prior = prRFreqMulti2;
          } else if (degree == 3) {
            prior = prRFreqMulti3;
          } else if (degree == 4) {
            prior = prRFreqMulti4;
          } else if (degree == 5) {
            prior = prRFreqMulti5;
          }
          break;
        case d_probit:
          printf("error with restricted probit model\n");
          return;
          break;
        case d_qlinear:
          printf("error with restricted QLinear model\n");
          return;
          break;
        case d_weibull:
          anal.model = d_weibull;
          prior = prRFreqWeibull;
          break;
        default:
          printf("error with restricted models\n");
          return;
      }
    } else {
      // unrestricted
      switch (model) {
        case d_hill:
          anal.model = d_hill;
          prior = prUFreqDHill;
          break;
        case d_gamma:
          anal.model = d_gamma;
          prior = prUFreqGamma;
          break;
        case d_logistic:
          anal.model = d_logistic;
          prior = prUFreqLogistic;
          break;
        case d_loglogistic:
          anal.model = d_loglogistic;
          prior = prUFreqLogLogistic;
          break;
        case d_logprobit:
          anal.model = d_logprobit;
          prior = prUFreqLogProbit;
          break;
        case d_multistage:
          anal.model = d_multistage;
          if (degree == 1) {
            prior = prUFreqMulti1;
          } else if (degree == 2) {
            prior = prUFreqMulti2;
          } else if (degree == 3) {
            prior = prUFreqMulti3;
          } else if (degree == 4) {
            prior = prUFreqMulti4;
          } else if (degree == 5) {
            prior = prUFreqMulti5;
          }
          break;
        case d_probit:
          anal.model = d_probit;
          prior = prUFreqProbit;
          break;
        case d_qlinear:
          anal.model = d_qlinear;
          prior = prUFreqQLinear;
          break;
        case d_weibull:
          anal.model = d_weibull;
          prior = prUFreqWeibull;
          break;
        default:
          printf("error with restricted models\n");
          return;
      }
    }
  } else {
    // bayesian
    switch (model) {
      case d_hill:
        anal.model = d_hill;
        prior = prBayesianDHill;
        break;
      case d_gamma:
        anal.model = d_gamma;
        prior = prBayesianGamma;
        break;
      case d_logistic:
        anal.model = d_logistic;
        prior = prBayesianLogistic;
        return;
        break;
      case d_loglogistic:
        anal.model = d_loglogistic;
        prior = prBayesianLogLogistic;
        break;
      case d_logprobit:
        anal.model = d_logprobit;
        prior = prBayesianLogProbit;
        break;
      case d_multistage:
        anal.model = d_multistage;
        if (degree == 1) {
          prior = prBayesianMulti1;
        } else if (degree == 2) {
          prior = prBayesianMulti2;
        } else if (degree == 3) {
          prior = prBayesianMulti3;
        } else if (degree == 4) {
          prior = prBayesianMulti4;
        } else if (degree == 5) {
          prior = prBayesianMulti5;
        }
        break;
      case d_probit:
        anal.model = d_probit;
        prior = prBayesianProbit;
        break;
      case d_qlinear:
        anal.model = d_qlinear;
        prior = prBayesianQLinear;
        break;
      case d_weibull:
        anal.model = d_weibull;
        prior = prBayesianWeibull;
        break;
      default:
        printf("error with restricted models\n");
        return;
    }
  }

  // declare analysis
  anal.BMD_type = BMD_type;
  anal.BMR = BMR;
  anal.alpha = alpha;
  anal.parms = numParms;
  anal.Y.assign(Y, Y + numDataRows);
  anal.n_group.assign(N, N + numDataRows);
  anal.doses.assign(D, D + numDataRows);
  anal.prior_cols = prCols;
  anal.n = numDataRows;
  anal.degree = degree;
  anal.prior.assign(prior, prior + anal.prior_cols * anal.parms);

  struct python_dichotomous_model_result res;
  res.model = anal.model;
  res.dist_numE = 200;
  res.nparms = anal.parms;

  pythonBMDSDicho(&anal, &res);

  std::cout << "INPUT" << std::endl;
  printBmdsStruct(&anal);

  std::cout << "OUTPUT" << std::endl;
  printBmdsStruct(&res);

  // printDichoModResult(&anal, &res, showResultsOverride);
}

void runDichoMA() {
  //
  //  printf("Running dichotomous Model Averaging\n");
  //
  ////estimate_ma_laplace_dicho(struct dichotomousMA_analysis *MA,
  ////                         struct dichotomous_analysis *DA ,
  ////                         struct dichotomousMA_result *res);
  //
  /////////////////////////////////
  ////USER INPUT
  /////////////////////////////////
  //
  //  #define numDichoModelsP1 10  //# of dicho models + 1
  //  #define numModels 9  //1-9
  //  bool include[numDichoModelsP1];
  //  include[d_hill]  = true;
  //  include[d_gamma] = true;
  //  include[d_logistic] = true;
  //  include[d_loglogistic] = true;
  //  include[d_logprobit] = true;
  //  include[d_multistage] = true;
  //  include[d_probit] = true;
  //  include[d_qlinear] = true;
  //  include[d_weibull] = true;
  //
  //
  //
  //  double modelPriors[numModels];
  //  modelPriors[0] = 1.0/numModels;
  //  modelPriors[1] = 1.0/numModels;
  //  modelPriors[2] = 1.0/numModels;
  //  modelPriors[3] = 1.0/numModels;
  //  modelPriors[4] = 1.0/numModels;
  //  modelPriors[5] = 1.0/numModels;
  //  modelPriors[6] = 1.0/numModels;
  //  modelPriors[7] = 1.0/numModels;
  //  modelPriors[8] = 1.0/numModels;
  ////  #define numModels 9
  ////  enum dich_model model = d_hill;  //d_hill =1, d_gamma=2,d_logistic=3, d_loglogistic=4,
  //                                   //d_logprobit=5, d_multistage=6,d_probit=7,
  //                                   //d_qlinear=8,d_weibull=9
  //  int BMD_type = 1;        // 1 = extra ; added otherwise
  //  double BMR = 0.1;
  //  double alpha = 0.05;
  /////////////////////////////////
  /////////////////////////////////
  ////dicho data - dose, N, incidence
  /////////////////////////////////
  //
  ////  double D[] = {0,50, 100, 150, 200};
  ////  double Y[] = {0, 5, 30, 65, 90};
  ////  double N[] = {100, 100, 100, 100, 100};
  //
  ////  user submitted
  //    double D[] = {0, 11, 30, 100, 356};
  //    double Y[] = {2, 10, 13, 15, 15};
  //    double N[] = {14, 15, 15, 15, 15};
  //
  ////  D80
  ////    double D[] = {0, 4.79, 9.57};
  ////    double Y[] = {0, 7, 8};
  ////    double N[] = {10, 15, 24};
  //
  ////  D100
  ////  double D[] = {0, 5.1, 21.9, 46.5};
  ////  double Y[] = {5, 5, 9, 17};
  ////  double N[] = {60, 60, 60, 60};
  //
  ////  D101
  ////  double D[] = {0, 1127, 2435, 5203};
  ////  double Y[] = {2, 8, 9, 30};
  ////  double N[] = {50, 49, 49, 50};
  //
  ////  D107
  ////  double D[] = {0, 209.8, 444.6, 978.1};
  ////  double Y[] = {8, 17, 26, 42};
  ////  double N[] = {50, 50, 50, 50};
  //
  ////  D110
  ////  double D[] = {0, 93.33, 196.4, 403.4};
  ////  double Y[] = {2, 2, 3, 8};
  ////  double N[] = {50, 50, 50, 50};
  //
  ////  D116
  ////  double D[] = {0, 0.156, 0.312, 0.625, 1.25, 2.5};
  ////  double Y[] = {0, 0, 0, 2, 10, 10};
  ////  double N[] = {10, 10, 10, 10, 10, 10};
  //
  //
  ///////////////////////////////////////////////////
  //////END USER INPUT
  //////////////////////////////////////////////////////
  //
  //  //check included models vs numModels
  //  int tmpNumModels = 0;
  //  for (int i=1; i<numDichoModelsP1; i++){
  //    if(include[i] == true){
  //      tmpNumModels++;
  //    }
  //  }
  //  if (tmpNumModels != numModels){
  //    printf("Number of models does not match included models\n");
  //    exit(-1);
  //  }
  //
  //  int numDataRows = sizeof(D)/sizeof(D[0]);
  //
  //  //check data array sizes for consistency
  //  size_t numElementsY = sizeof(Y)/sizeof(Y[0]);
  //  size_t numElementsN = sizeof(N)/sizeof(N[0]);
  //  if (numDataRows != numElementsY || numElementsY != numElementsN) {
  //    printf("Number of data elements are not consistent\nExiting Code\n");
  //    exit(-1);
  //  }
  //
  //  //create array of priors
  //
  //  struct dichotomous_analysis anal;
  //  anal.BMD_type = BMD_type;
  //  anal.BMR = BMR;
  //  anal.alpha = alpha;
  //  anal.Y = Y;
  //  anal.n_group = N;
  //  anal.doses = D;
  //  anal.n = numDataRows;
  ////  anal.samples = ????
  ////  anal.burnin = ?????
  //
  //
  //  printf("numModels: %d\n", numModels);
  //
  //  //priors defined columnwise
  //  int prCols = 5;
  //
  //  //run selected models in MA
  ////  int models[numModels];
  ////  int priorCols[numModels];
  ////  for (int i=0; i<numModels; i++){
  //////    models[i] = i+1;
  ////    priorCols[i] = prCols;
  ////  }
  //
  //
  //  //int numParms[numModels];
  //
  //  std::vector<int> models(numModels);
  //  std::vector<int> priorCols(numModels, prCols);
  //  std::vector<int> numParms(numModels);
  //  //double prHill[20] = {1,1,1,2,-1,0,-3,0.693147,2,3,3.3,0.5,-40,-40,-40,0,40,40,40,40};
  //  //double prGamma[15] = {1,2,2,0,0.693147,0,2,0.424264,1,-18,0.2,0,18,20,10000};
  //
  //  //double prLog[10] = {1,2,0,0,2,2,-20,0,20,40};
  //  //double prLogLog[15] = {1,1,2,0,0,0.693147,2,1,0.5,-20,-40,0,20,40,20}; //d_loglogistic
  //  //double prLogProb[15] = {1,1,2,0,0,0.693147,2,1,0.5,-20,-40,0,20,40,20}; //d_logprobit
  //  //double prMulti[20] = {1,2,2,2,0,0,0,0,2,1,1,1,-20,0,0,0,20,1e6,1e6,1e6}; //d_multistage
  //  //double prProb[10] = {1,2,0,0,1,2,-20,0,20,40}; //d_probit
  //  //double prQLin[10] = {1,2,0,0,2,1,-20,0,20,18}; //d_qlinear
  //  //double prWeib[15] = {1,2,2,0,0.424264,0,2,0.5,1.5,-20,0,0,20,40,1e5}; //d_weibull
  ////  double prHill[20] = {1,1,1,2,-1,0,-3,0.693147,2,3,3.3,0.5,-40,-40,-40,1e-8,40,40,40,40};
  ////  double prGamma[15] = {1,2,2,0,0.693147,0,2,0.424264,1,-18,0.2,1e-4,18,20,100};
  ////
  ////  double prLog[10] = {1,2,0,0.1,2,1,-20,1e-12,20,100};
  ////  double prLogLog[15] = {1,1,2,0,0,0.693147,2,1,0.5,-20,-40,1e-4,20,40,20}; //d_loglogistic
  ////  double prLogProb[15] = {1,1,2,0,0,0.693147,2,1,0.5,-20,-8,1e-4,20,8,40}; //d_logprobit
  ////  double prMulti[20] = {1,2,2,2,0,0,0,0,2,0.5,1,1,-20,1e-4,1e-4,1e-4,20,100,1e6,1e6};
  /////d_multistage /  double prProb[10] = {1,2,0,0.1,2,1,-8,0,8,40}; //d_probit /  double
  /// prQLin[10] = {1,2,0,0.5,2,1,-20,0,20,100}; //d_qlinear /  double prWeib[15] =
  ///{1,2,2,0,0.693147,0,2,0.424264,1,-20,1e-4,1e-4,20,18,20}; //d_weibull
  //
  //  //double *pr[numModels] = {prHill,prGamma};
  //  double *pr[numModels];
  //  //pr[0] = prHill;
  //  //pr[1] = prGamma;
  //
  //  int count = 0;
  //  if(include[d_hill]){
  //    //pr[count] = prHill;
  //    //numParms[count] = sizeof(prHill)/sizeof(prHill[0])/prCols;
  //    pr[count] = prBayesianDHill;
  //    numParms[count] = sizeof(prBayesianDHill)/sizeof(prBayesianDHill[0])/prCols;
  //    models[count] = d_hill;
  //    count++;
  //  }
  //  if(include[d_gamma]){
  //    //pr[count] = prGamma;
  //    //numParms[count] = sizeof(prGamma)/sizeof(prGamma[0])/prCols;
  //    pr[count] = prBayesianGamma;
  //    numParms[count] = sizeof(prBayesianGamma)/sizeof(prBayesianGamma[0])/prCols;
  //    models[count] = d_gamma;
  //    count++;
  //  }
  //  if(include[d_logistic]){
  //    //pr[count] = prLog;
  //    //numParms[count] = sizeof(prLog)/sizeof(prLog[0])/prCols;
  //    pr[count] = prBayesianLogistic;
  //    numParms[count] = sizeof(prBayesianLogistic)/sizeof(prBayesianLogistic[0])/prCols;
  //    models[count] = d_logistic;
  //    count++;
  //  }
  //  if(include[d_loglogistic]){
  //    //pr[count] = prLogLog;
  //    //numParms[count] = sizeof(prLogLog)/sizeof(prLogLog[0])/prCols;
  //    pr[count] = prBayesianLogLogistic;
  //    numParms[count] = sizeof(prBayesianLogLogistic)/sizeof(prBayesianLogLogistic[0])/prCols;
  //    models[count] = d_loglogistic;
  //    count++;
  //  }
  //  if(include[d_logprobit]){
  //    //pr[count] = prLogProb;
  //    //numParms[count] = sizeof(prLogProb)/sizeof(prLogProb[0])/prCols;
  //    pr[count] = prBayesianLogProbit;
  //    numParms[count] = sizeof(prBayesianLogProbit)/sizeof(prBayesianLogProbit[0])/prCols;
  //    models[count] = d_logprobit;
  //    count++;
  //  }
  //  if(include[d_multistage]){
  //    //pr[count] = prMulti;
  //    //numParms[count] = sizeof(prMulti)/sizeof(prMulti[0])/prCols;
  //    pr[count] = prBayesianMulti3;
  //    numParms[count] = sizeof(prBayesianMulti3)/sizeof(prBayesianMulti3[0])/prCols;
  //    models[count] = d_multistage;
  //    count++;
  //  }
  //  if(include[d_probit]){
  //    //pr[count] = prProb;
  //    //numParms[count] = sizeof(prProb)/sizeof(prProb[0])/prCols;
  //    pr[count] = prBayesianProbit;
  //    numParms[count] = sizeof(prBayesianProbit)/sizeof(prBayesianProbit[0])/prCols;
  //    models[count] = d_probit;
  //    count++;
  //  }
  //  if(include[d_qlinear]){
  //    //pr[count] = prQLin;
  //    //numParms[count] = sizeof(prQLin)/sizeof(prQLin[0])/prCols;
  //    pr[count] = prBayesianQLinear;
  //    numParms[count] = sizeof(prBayesianQLinear)/sizeof(prBayesianQLinear[0])/prCols;
  //    models[count] = d_qlinear;
  //    count++;
  //  }
  //  if(include[d_weibull]){
  //    //pr[count] = prWeib;
  //    //numParms[count] = sizeof(prWeib)/sizeof(prWeib[0])/prCols;
  //    pr[count] = prBayesianWeibull;
  //    numParms[count] = sizeof(prBayesianWeibull)/sizeof(prBayesianWeibull[0])/prCols;
  //    models[count] = d_weibull;
  //    count++;
  //  }
  //
  ////  int count = 0;
  ////  numParms[count] = sizeof(prHill)/sizeof(prHill[0])/prCols;
  ////  models[count] = 1;
  ////  count++;
  ////  numParms[count] = sizeof(prGamma)/sizeof(prGamma[0])/prCols;
  ////  models[count] = 2;
  ////  count++;
  ////  numParms[count] = sizeof(prLog)/sizeof(prLog[0])/prCols;
  ////  count++;
  ////  numParms[count] = sizeof(prLogLog)/sizeof(prLogLog[0])/prCols;
  ////  count++;
  ////  numParms[count] = sizeof(prLogProb)/sizeof(prLogProb[0])/prCols;
  ////  count++;
  ////  numParms[count] = sizeof(prMulti)/sizeof(prMulti[0])/prCols;
  ////  count++;
  ////  numParms[count] = sizeof(prProb)/sizeof(prProb[0])/prCols;
  ////  count++;
  ////  numParms[count] = sizeof(prQLin)/sizeof(prQLin[0])/prCols;
  ////  count++;
  ////  numParms[count] = sizeof(prWeib)/sizeof(prWeib[0])/prCols;
  ////  count++;
  //
  //  if (count != numModels) {
  //    printf("Error in specifying parameters");
  //    return;
  //  }
  //
  //  for (int i=0; i<numModels; i++){
  //    printf("model %d has %d parms\n",i,numParms[i]);
  //  }
  //
  //
  //  //parms array declared
  //
  ////  printf("numModels=%d\n",numModels);
  ////  printf("prCols=%d\n",prCols);
  //  printf("model priors1\n");
  //  for (int i=0; i<numModels; i++){
  //    printf("Model:%d\n",i);
  //    for (int j=0; j<numParms[i]*prCols; j++){
  //      printf("%f, ", *(pr[i] + j));
  //    }
  //    printf("\n");
  //  }
  //
  //
  //  struct python_dichotomousMA_analysis ma_info;
  //  ma_info.actual_parms = numParms;
  //  ma_info.prior_cols = priorCols;
  //  ma_info.models = models;
  //  ma_info.priors = pr;
  //  ma_info.modelPriors = modelPriors;
  //  ma_info.nmodels = numModels;
  //
  //  struct dichotomous_model_result *res[numModels];
  //  int dist_numE = 200;
  //
  //
  //  //for (int i=0; i<numModels; i++){
  //  //  struct dichotomous_model_result modelRes;
  //  //  modelRes.model = models[i];
  //  //  modelRes.nparms = numParms[i];
  //  //  modelRes.dist_numE = dist_numE;
  //  //  modelRes.parms = malloc(sizeof(double)*numParms[i]);
  //  //  modelRes.cov = malloc(sizeof(double)*numParms[i]*numParms[i]);
  //  //  modelRes.bmd_dist = malloc(sizeof(double)*dist_numE*2);
  //  //  res[i] = &modelRes;
  //  //}
  //
  //
  //  for (int i=0; i<numModels; i++){
  ////    res[i] = malloc(sizeof(struct dichotomous_model_result));
  //    res[i] = new dichotomous_model_result;
  //    res[i]->model = models[i];
  //    res[i]->nparms = numParms[i];
  //    res[i]->dist_numE = dist_numE;
  //    res[i]->parms = (double*)malloc(sizeof(double)*numParms[i]);
  //    res[i]->cov = (double*)malloc(sizeof(double)*numParms[i]*numParms[i]);
  //    res[i]->bmd_dist = (double*)malloc(sizeof(double)*dist_numE*2);
  //    //struct dichotomous_model_result modelRes;
  //    //modelRes.model = models[i];
  //    //modelRes.nparms = numParms[i];
  //    //modelRes.dist_numE = dist_numE;
  //    //modelRes.parms = malloc(sizeof(double)*numParms[i]);
  //    //modelRes.cov = malloc(sizeof(double)*numParms[i]*numParms[i]);
  //    //modelRes.bmd_dist = malloc(sizeof(double)*dist_numE*2);
  //    //res[i] = &modelRes;
  //  }
  //
  ////    struct dichotomous_model_result modelRes0;
  ////    modelRes0.model = models[0];
  ////    modelRes0.nparms = numParms[0];
  ////    modelRes0.dist_numE = dist_numE;
  ////    modelRes0.parms = malloc(sizeof(double)*numParms[0]);
  ////    modelRes0.cov = malloc(sizeof(double)*numParms[0]*numParms[0]);
  ////    modelRes0.bmd_dist = malloc(sizeof(double)*dist_numE*2);
  ////    res[0] = &modelRes0;
  ////    struct dichotomous_model_result modelRes1;
  ////    modelRes1.model = models[1];
  ////    modelRes1.nparms = numParms[1];
  ////    modelRes1.dist_numE = dist_numE;
  ////    modelRes1.parms = malloc(sizeof(double)*numParms[1]);
  ////    modelRes1.cov = malloc(sizeof(double)*numParms[1]*numParms[1]);
  ////    modelRes1.bmd_dist = malloc(sizeof(double)*dist_numE*2);
  ////    res[1] = &modelRes1;
  //  struct dichotomousMA_result ma_res;
  //  ma_res.nmodels = numModels;
  //  ma_res.models = res;
  //  ma_res.dist_numE = dist_numE;
  //  //double post_probs[numModels];
  //  double* post_probs = new double[numModels];
  //  ma_res.post_probs = post_probs;
  //  //double bmd_dist[dist_numE*2];
  //  double* bmd_dist = new double[dist_numE*2];
  //   ma_res.bmd_dist = bmd_dist;
  //
  //
  ////  for (int i=0; i<numModels; i++){
  ////    printf("dist_numE=%d\n",ma_res.models[i]->dist_numE);
  ////  }
  //
  ////  printf("calling estimate_ma_laplace\n");
  ////  estimate_ma_laplace_dicho(&ma_info, &anal, &ma_res);
  //
  //
  //  struct BMDSMA_results bmdsRes;
  //
  //  //double BMD[numModels];
  //  //double BMDL[numModels];
  //  //double BMDU[numModels];
  //  //double ebLower[anal.n];
  //  //double ebUpper[anal.n];
  //  //double* BMD = new double[numModels];
  //  //double* BMDL = new double[numModels];
  //  //double* BMDU = new double[numModels];
  //  //double* ebLower = new double[anal.n];
  //  //double* ebUpper = new double[anal.n];
  //  std::vector<double> BMD(numModels, BMDS_MISSING);
  //  std::vector<double> BMDL(numModels, BMDS_MISSING);
  //  std::vector<double> BMDU(numModels, BMDS_MISSING);
  //  std::vector<double> ebLower(anal.n, BMDS_MISSING);
  //  std::vector<double> ebUpper(anal.n, BMDS_MISSING);
  //
  //  //for (int i=0; i<numModels; i++){
  //  //  BMD[i] = -9999.0;
  //  //  BMDL[i] = -9999.0;
  //  //  BMDU[i] = -9999.0;
  //  //}
  //  //for (int i=0; i<anal.n; i++){
  //  //  ebLower[i] = -9999.0;
  //  //  ebUpper[i] = -9999.0;
  //  //}
  //  bmdsRes.BMD_MA = BMDS_MISSING;
  //  bmdsRes.BMDU_MA = BMDS_MISSING;
  //  bmdsRes.BMDL_MA = BMDS_MISSING;
  //  bmdsRes.BMD = BMD;
  //  bmdsRes.BMDL = BMDL;
  //  bmdsRes.BMDU = BMDU;
  //  bmdsRes.ebLower = ebLower;
  //  bmdsRes.ebUpper = ebUpper;
  //
  //  //runBMDSDichoMA(&ma_info, &anal, &ma_res, &bmdsRes);
  //  pythonBMDSDichoMA(&ma_info, &anal, &ma_res, &bmdsRes);
  //
  //
  ////  printf("\nBMD Dist:\n");
  ////  for (int i=0; i<ma_res.dist_numE; i++){
  ////    printf("i:%d, perc:%f, dist:%f\n", i, ma_res.bmd_dist[i+ma_res.dist_numE],
  /// ma_res.bmd_dist[i]); /  } /  struct dichotomous_model_result indRes;
  //
  ////  printf("individual BMD Dist:\n");
  ////  for (int j=0; j<numModels; j++){
  ////    indRes = *ma_res.models[j];
  ////    printf("\nModel %d\n", j);
  ////    for (int i=0; i<indRes.dist_numE; i++){
  ////      printf("i:%d, perc:%f, dist:%f\n", i, indRes.bmd_dist[i+indRes.dist_numE],
  /// indRes.bmd_dist[i]); /    } /  }
  //
  //
  //  printf("\nBenchmark Dose\n");
  //  printf("MA BMD: %f\n",bmdsRes.BMD_MA);
  //  printf("MA BMDL: %f\n",bmdsRes.BMDL_MA);
  //  printf("MA BMDU: %f\n",bmdsRes.BMDU_MA);
  //
  //  printf("\nMA - Individual Models\n");
  //  for(int i=0; i<numModels; i++){
  //    printf("i:%d, model:%d\n", i, ma_res.models[i]->model);
  //    printf("\tpost prob:%f\n", ma_res.post_probs[i]);
  //    printf("\tBMD:%f\n",bmdsRes.BMD[i]);
  //    printf("\tBMDL:%f\n",bmdsRes.BMDL[i]);
  //    printf("\tBMDU:%f\n",bmdsRes.BMDU[i]);
  //    printf("\tParms:\n");
  //    for(int j=0; j<ma_res.models[i]->nparms; j++){
  //      printf("\t\tj:%d, value:%f\n", j, ma_res.models[i]->parms[j]);
  //    }
  //    //printf("i:%d, model:%d, post prob:%f, BMD:%f, BMDL:%f,
  //    BMDU:%f\n",i,ma_res.models[i]->model,ma_res.post_probs[i],bmdsRes.BMD[i],bmdsRes.BMDL[i],bmdsRes.BMDU[i]);
  //  }
  //  printf("Error bars\n");
  //  for(int i=0; i<anal.n; i++){
  //    printf("%f\t%f\n", ebLower[i], ebUpper[i]);
  //  }
  //
  ////  printf("\nMA - Individual Models\n");
  ////  for(int i=0; i<numModels; i++){
  ////    printf("i:%d, model:%d, post prob:%f, BMD:%f, BMDL:%f,
  /// BMDU:%f\n",i,ma_res.models[i]->model,ma_res.post_probs[i],bmdsRes.BMD[i],bmdsRes.BMDL[i],bmdsRes.BMDU[i]);
  ////  }
  //
  ////  printf("model priors2\n");
  ////  for (int i=0; i<numModels; i++){
  ////    printf("model %d\n", i);
  ////    for (int j=0; j<numParms[i]*prCols; j++){
  ////      printf("%f, ", *(ma_info.priors[i] + j));
  ////    }
  ////    printf("\n");
  ////  }
}

void runPythonDichoMA() {
  printf("Running dichotomous Model Averaging\n");

  // estimate_ma_laplace_dicho(struct dichotomousMA_analysis *MA,
  //                          struct dichotomous_analysis *DA ,
  //                          struct dichotomousMA_result *res);

  ///////////////////////////////
  // USER INPUT
  ///////////////////////////////

#define numDichoModelsP1 10  // # of dicho models + 1
#define numModels 9          // 1-9
  //  bool include[numDichoModelsP1];
  //  include[d_hill]  = true;
  //  include[d_gamma] = true;
  //  include[d_logistic] = true;
  //  include[d_loglogistic] = true;
  //  include[d_logprobit] = true;
  //  include[d_multistage] = true;
  //  include[d_probit] = true;
  //  include[d_qlinear] = true;
  //  include[d_weibull] = true;
  std::vector<bool> include(numModels, true);

  //  double modelPriors[numModels];
  //  modelPriors[0] = 1.0/numModels;
  //  modelPriors[1] = 1.0/numModels;
  //  modelPriors[2] = 1.0/numModels;
  //  modelPriors[3] = 1.0/numModels;
  //  modelPriors[4] = 1.0/numModels;
  //  modelPriors[5] = 1.0/numModels;
  //  modelPriors[6] = 1.0/numModels;
  //  modelPriors[7] = 1.0/numModels;
  //  modelPriors[8] = 1.0/numModels;

  std::vector<double> modelPriors(numModels, 1.0 / numModels);
  //  #define numModels 9
  //  enum dich_model model = d_hill;  //d_hill =1, d_gamma=2,d_logistic=3, d_loglogistic=4,
  // d_logprobit=5, d_multistage=6,d_probit=7,
  // d_qlinear=8,d_weibull=9
  int BMD_type = 1;  // 1 = extra ; added otherwise
  double BMR = 0.1;
  double alpha = 0.05;
  ///////////////////////////////
  ///////////////////////////////
  // dicho data - dose, N, incidence
  ///////////////////////////////

  //  double D[] = {0,50, 100, 150, 200};
  //  double Y[] = {0, 5, 30, 65, 90};
  //  double N[] = {100, 100, 100, 100, 100};

  //  user submitted
  double D[] = {0, 11, 30, 100, 356};
  double Y[] = {2, 10, 13, 15, 15};
  double N[] = {14, 15, 15, 15, 15};

  //  D80
  //    double D[] = {0, 4.79, 9.57};
  //    double Y[] = {0, 7, 8};
  //    double N[] = {10, 15, 24};

  //  D100
  //  double D[] = {0, 5.1, 21.9, 46.5};
  //  double Y[] = {5, 5, 9, 17};
  //  double N[] = {60, 60, 60, 60};

  //  D101
  //  double D[] = {0, 1127, 2435, 5203};
  //  double Y[] = {2, 8, 9, 30};
  //  double N[] = {50, 49, 49, 50};

  //  D107
  //  double D[] = {0, 209.8, 444.6, 978.1};
  //  double Y[] = {8, 17, 26, 42};
  //  double N[] = {50, 50, 50, 50};

  //  D110
  //  double D[] = {0, 93.33, 196.4, 403.4};
  //  double Y[] = {2, 2, 3, 8};
  //  double N[] = {50, 50, 50, 50};

  //  D116
  //  double D[] = {0, 0.156, 0.312, 0.625, 1.25, 2.5};
  //  double Y[] = {0, 0, 0, 2, 10, 10};
  //  double N[] = {10, 10, 10, 10, 10, 10};

  /////////////////////////////////////////////////
  ////END USER INPUT
  ////////////////////////////////////////////////////

  // check included models vs numModels
  int tmpNumModels = 0;
  for (int i = 1; i < numDichoModelsP1; i++) {
    if (include[i] == true) {
      tmpNumModels++;
    }
  }
  if (tmpNumModels != numModels) {
    printf("Number of models does not match included models\n");
    exit(-1);
  }

  int numDataRows = sizeof(D) / sizeof(D[0]);

  // check data array sizes for consistency
  size_t numElementsY = sizeof(Y) / sizeof(Y[0]);
  size_t numElementsN = sizeof(N) / sizeof(N[0]);
  if (numDataRows != numElementsY || numElementsY != numElementsN) {
    printf("Number of data elements are not consistent\nExiting Code\n");
    exit(-1);
  }

  // create array of priors

  struct python_dichotomous_analysis anal;
  anal.BMD_type = BMD_type;
  anal.BMR = BMR;
  anal.alpha = alpha;
  anal.Y.assign(Y, Y + numDataRows);
  anal.n_group.assign(N, N + numDataRows);
  anal.doses.assign(D, D + numDataRows);
  anal.n = numDataRows;

  printf("numModels: %d\n", numModels);

  // priors defined columnwise
  int prCols = 5;

  std::vector<int> models;
  std::vector<int> priorCols(numModels, prCols);
  std::vector<int> numParms;
  std::vector<std::vector<double>> pr;
  double *prArray;
  std::vector<double> curPR;
  int prSize;

  //  int count = 0;
  if (include[d_hill]) {
    models.push_back(d_hill);
    prArray = prBayesianDHill;
    prSize = sizeof(prBayesianDHill) / sizeof(prBayesianDHill[0]);
    curPR.assign(prArray, prArray + prSize);
    pr.push_back(curPR);
    numParms.push_back(prSize / prCols);
  }
  if (include[d_gamma]) {
    models.push_back(d_gamma);
    prArray = prBayesianGamma;
    prSize = sizeof(prBayesianGamma) / sizeof(prBayesianGamma[0]);
    curPR.assign(prArray, prArray + prSize);
    pr.push_back(curPR);
    numParms.push_back(prSize / prCols);
  }
  if (include[d_logistic]) {
    models.push_back(d_logistic);
    prArray = prBayesianLogistic;
    prSize = sizeof(prBayesianLogistic) / sizeof(prBayesianLogistic[0]);
    curPR.assign(prArray, prArray + prSize);
    pr.push_back(curPR);
    numParms.push_back(prSize / prCols);
  }
  if (include[d_loglogistic]) {
    models.push_back(d_loglogistic);
    prArray = prBayesianLogLogistic;
    prSize = sizeof(prBayesianLogLogistic) / sizeof(prBayesianLogLogistic[0]);
    curPR.assign(prArray, prArray + prSize);
    pr.push_back(curPR);
    numParms.push_back(prSize / prCols);
  }
  if (include[d_logprobit]) {
    models.push_back(d_logprobit);
    prArray = prBayesianLogProbit;
    prSize = sizeof(prBayesianLogProbit) / sizeof(prBayesianLogProbit[0]);
    curPR.assign(prArray, prArray + prSize);
    pr.push_back(curPR);
    numParms.push_back(prSize / prCols);
  }
  if (include[d_multistage]) {
    models.push_back(d_multistage);
    prArray = prBayesianMulti3;
    prSize = sizeof(prBayesianMulti3) / sizeof(prBayesianMulti3[0]);
    curPR.assign(prArray, prArray + prSize);
    pr.push_back(curPR);
    numParms.push_back(prSize / prCols);
  }
  if (include[d_probit]) {
    models.push_back(d_probit);
    prArray = prBayesianProbit;
    prSize = sizeof(prBayesianProbit) / sizeof(prBayesianProbit[0]);
    curPR.assign(prArray, prArray + prSize);
    pr.push_back(curPR);
    numParms.push_back(prSize / prCols);
  }
  if (include[d_qlinear]) {
    models.push_back(d_qlinear);
    prArray = prBayesianQLinear;
    prSize = sizeof(prBayesianQLinear) / sizeof(prBayesianQLinear[0]);
    curPR.assign(prArray, prArray + prSize);
    pr.push_back(curPR);
    numParms.push_back(prSize / prCols);
  }
  if (include[d_weibull]) {
    models.push_back(d_weibull);
    prArray = prBayesianWeibull;
    prSize = sizeof(prBayesianWeibull) / sizeof(prBayesianWeibull[0]);
    curPR.assign(prArray, prArray + prSize);
    pr.push_back(curPR);
    numParms.push_back(prSize / prCols);
  }

  if (models.size() != numModels) {
    printf("Error in specifying parameters");
    return;
  }

  for (int i = 0; i < numModels; i++) {
    printf("model %d has %d parms\n", i, numParms[i]);
  }

  struct python_dichotomousMA_analysis ma_info;
  ma_info.actual_parms = numParms;
  ma_info.prior_cols = priorCols;
  ma_info.models = models;
  ma_info.priors = pr;
  ma_info.modelPriors = modelPriors;
  ma_info.nmodels = numModels;
  ma_info.pyDA = anal;

  std::vector<python_dichotomous_model_result> res(numModels);
  int dist_numE = 200;

  for (int i = 0; i < numModels; i++) {
    res[i].model = models[i];
    res[i].nparms = numParms[i];
    res[i].dist_numE = dist_numE;
  }

  struct python_dichotomousMA_result ma_res;
  ma_res.nmodels = numModels;
  ma_res.models = res;
  ma_res.dist_numE = dist_numE;

  struct BMDSMA_results bmdsRes;
  bmdsRes.BMD.assign(numModels, BMDS_MISSING);
  bmdsRes.BMDL.assign(numModels, BMDS_MISSING);
  bmdsRes.BMDU.assign(numModels, BMDS_MISSING);
  bmdsRes.ebLower.assign(anal.n, BMDS_MISSING);
  bmdsRes.ebUpper.assign(anal.n, BMDS_MISSING);
  bmdsRes.BMD_MA = BMDS_MISSING;
  bmdsRes.BMDL_MA = BMDS_MISSING;
  bmdsRes.BMDU_MA = BMDS_MISSING;

  ma_res.bmdsRes = bmdsRes;

  pythonBMDSDichoMA(&ma_info, &ma_res);

  std::cout << "INPUT" << std::endl;
  printBmdsStruct(&ma_info);
  std::cout << "OUTPUT" << std::endl;
  printBmdsStruct(&ma_res);
}

void runPythonContLoud() {
  printf("Running LOUD Continuous Model Averaging\n");

  double D[] = {0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
                0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.125, 0.125,
                0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
                0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.250, 0.250, 0.250, 0.250,
                0.250, 0.250, 0.250, 0.250, 0.250, 0.250, 0.250, 0.250, 0.250, 0.250, 0.250,
                0.250, 0.250, 0.250, 0.250, 0.250, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500,
                0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500,
                0.500, 0.500, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
                1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000};

  double Yarr[] = {10.714175, 10.604914, 11.499425, 10.370608, 10.986296, 9.946773,  11.385419,
                   8.927492,  11.621717, 11.316601, 11.370181, 8.840015,  12.345874, 10.887448,
                   11.103332, 10.044302, 9.772624,  8.607546,  11.913568, 12.526342, 11.519071,
                   11.413047, 11.592883, 11.479883, 11.423585, 11.114481, 11.211780, 10.366999,
                   12.460339, 9.809126,  13.478086, 13.691934, 12.140507, 10.529764, 11.727003,
                   11.606124, 10.454382, 11.606850, 9.121244,  10.787093, 12.089929, 12.618154,
                   13.106158, 11.351711, 13.694759, 13.884344, 14.477708, 13.440013, 9.501834,
                   12.501327, 11.813062, 13.369222, 12.868868, 13.576607, 13.850305, 10.711790,
                   11.874109, 12.463099, 12.720392, 13.962033, 15.489551, 14.158361, 13.430108,
                   13.406044, 12.705685, 14.834740, 13.650278, 14.131006, 16.240119, 14.425319,
                   14.139144, 11.824363, 15.057471, 14.657698, 12.203581, 14.119994, 13.462503,
                   16.076048, 15.963720, 15.164853, 14.092978, 16.670621, 15.308121, 16.313754,
                   16.966145, 16.227856, 13.984670, 15.997495, 15.537613, 12.413045, 15.466326,
                   14.302947, 16.061184, 16.060429, 13.971971, 18.086398, 15.485381, 16.729798};

  bool suffStat = false;
  // bool isIncreasing = true;
  bool detectAdvDir = true;

  int numDataRows = sizeof(D) / sizeof(D[0]);

  // check data array sizes for consistency
  size_t numElementsY = sizeof(Yarr) / sizeof(Yarr[0]);
  if (numDataRows != numElementsY) {
    printf("Number of data elements are not consistent\nExiting Code\n");
    exit(-1);
  }

  std::vector<double> doses(D, D + numDataRows);
  std::vector<double> Y(Yarr, Yarr + numDataRows);

  // declare analysis
  struct python_continuous_analysis anal;
  anal.doses = doses;
  anal.Y = Y;
  anal.n = numDataRows;
  anal.detectAdvDir = detectAdvDir;

  //  if (suffStat) {
  //    anal.n_group.assign(N, N + numDataRows);
  //    anal.sd.assign(SD, SD + numDataRows);
  //  }
  //  anal.disttype = dist;
  //  if (!detectAdvDir) {
  //    anal.isIncreasing = isIncreasing;
  //  }

  //  anal.alpha = alpha;
  //  anal.BMD_type = BMD_type;  // 1=absdev, 2 = stddev, 3 = reldev, 4 = pt, 5 = extra, 6 =
  // hybrid_extra, 7 = hybrid_added   from src/include/cmodeldefs.h
  //  anal.BMR = BMRF;
  //  anal.samples = 0;  // num MCMC samples
  //  anal.tail_prob = 0.01;
  anal.suff_stat = suffStat;
  //  anal.isIncreasing = isIncreasing;
  //  anal.parms = numParms;
  //  anal.prior_cols = prCols;
  //  anal.transform_dose = 0;
  //  anal.prior.assign(prior, prior + anal.prior_cols * anal.parms);
  //  anal.restricted = restricted;
  //  anal.detectAdvDir = detectAdvDir;

  struct python_continuousMA_analysis ma_info;
  // ma_info.actual_parms = numParms;
  // ma_info.prior_cols = priorCols;
  // ma_info.models = models;
  // ma_info.priors = pr;
  // ma_info.modelPriors = modelPriors;
  // ma_info.nmodels = numModels;
  ma_info.pyCA = anal;

  struct python_continuousMA_result ma_res;

  pythonBMDSContLoud(&ma_info, &ma_res);
}

void runOldContAnalysis() {
  printf("Running continuous analysis\n");

  // char * bmdsVersion[32];

  // version(bmdsVersion);
  string bmdsVersion = version();
  // printf("Version: %s\n", bmdsVersion);
  std::cout << "Version: " << bmdsVersion << std::endl;

  bool isIncreasing;

  ///////////////////////////////
  // USER INPUT
  //////////////////////////////

  enum cont_model model = exp_5;    // hill, exp_3, exp_5, power, funl, polynomial
  int modelType = 1;                // 1 = frequentist, 2 = bayesian
  bool restricted = true;           // only used for frequentist models
  enum distribution dist = normal;  // normal, normal_ncv, log_normal
  bool detectAdvDir = true;         // if false then need to set isIncreasing
  bool countAllParmsOnBoundary = true;
  // isIncreasing = true;

  int degree = 2;  // for polynomial only

  double alpha = 0.05;
  double BMRF = 1.0;  // 1.0;
  int BMD_type = 2;   // 1=absdev, 2 = stddev, 3 = reldev, 4 = pt, 5 = extra, 6 = hybrid_extra, 7 =
                      // hybrid_added   from src/include/cmodeldefs.h
  ////////////////////////////////////////////////
  // cont data - suff stat: dose, Y, N, SD
  // cont data - individual: dose, response
  /////////////////////////////////////////////
  bool suffStat = true;

  // continuous1.dax
  //  double D[] = {0,25,50, 100, 200};
  //  double Y[] = {6.0, 5.2, 2.4, 1.1, 0.75};
  //  double N[] = {20, 20, 19, 20, 20};
  //  double SD[] = {1.2, 1.1, 0.81, 0.74, 0.66};
  // isIncreasing = false;

  // continuous2.dax
  //  double D[] = {0,0,0,0,18,18,18,18,18,20,20,20,20,30,30,30,30,35,35,35,35,40,40,40,40,40};
  //  double Y[] =
  //  {39,38.4,36.3,37.1,40.2,45.3,42.1,38.3,35.9,42.5,45.2,40.1,39.8,50.1,53.4,48.2,52.1,56.1,50.4,53.2,55.2,55.1,59.1,56.3,52.9,53.7};
  //  double N[1];
  //  double SD[1];
  //  isIncreasing = true;

  // continuous3.dax
  //  double D[] = {0,35,105,316,625};
  //  double Y[] = {1.61,1.66,1.75,1.81,1.89};
  //  double N[] = {10,10,10,10,10};
  //  double SD[] = {0.12,0.13,0.11,0.15,0.13};

  // other test datasets
  //  double D[] = {0,50, 100, 150, 200};
  //  double Y[] = {10, 20 , 30, 40 ,50};
  //  double N[] = {100, 100, 100, 100, 100};
  //  double SD[] = {3, 4, 5, 6, 7};
  //  isIncreasing = true;

  //  double D[] = {0,50, 100, 150, 200};
  //  double Y[] = {10, 0 , -10, -20 ,-30};
  //  double N[] = {100, 100, 100, 100, 100};
  //  double SD[] = {3, 4, 5, 6, 7};
  //  isIncreasing = false;

  //  double D[] = {0,50, 100, 150, 200};
  //  double Y[] = {10, 18, 32, 38, 70};
  //  double N[] = {100, 100, 100, 100, 100};
  //  double SD[] = {3.2, 4.8, 6.5, 7.2, 8.4};
  //  isIncreasing = true;

  //  double D[] = {0,50, 100, 150, 200};
  //  double Y[] = {1, -5 , -10, -20 ,-30};
  //  double N[] = {100, 100, 100, 100, 100};
  //  double SD[] = {3, 4, 5, 6, 7};
  //  isIncreasing = false;

  // c1b
  //  double D[] = {0, 75, 250};
  //  double Y[] = {15.81, 17.91, 21.48};
  //  double N[] = {11, 11, 11};
  //  double SD[] = {2.793, 2.902, 5.771};
  //  isIncreasing = true;

  // c2
  //    double D[] = {0,75,250};
  //    double Y[] = {94.18, 91.07, 116.61};
  //    double N[] = {11,11,11};
  //    double SD[] = {14.53, 9.22, 18.31};

  // c10
  //  double D[] = {0, 75, 250};
  //  double Y[] = {3.8, 8.9, 8.9};
  //  double N[] = {22, 22, 22};
  //  double SD[] = {3.34, 5.16, 2.49};

  // c20
  //  double D[] = {0, 1, 3, 9};
  //  double Y[] = {1.037, 1.05, 1.052, 1.066};
  //  double N[] = {10, 10, 10, 10};
  //  double SD[] = {0.015, 0.01, 0.01, 0.01};

  // c40
  //  double D[] = {0, 25, 100, 400};
  //  double Y[] = {0.67, 0.68, 0.71, 0.62};
  //  double N[] = {14, 15, 15, 15};
  //  double SD[] = {0.13, 0.09, 0.11, 0.08};

  // c60
  //  double D[] = {0, 10, 50, 100, 250};
  //  double Y[] = {0.116, 0.113, 0.108, 0.108, 0.106};
  //  double N[] = {30, 30, 30, 30, 30};
  //  double SD[] = {0.006, 0.006, 0.004, 0.009, 0.008};

  // c70b
  //  double D[] = {0, 46.4, 68.1, 200};
  //  double Y[] = {6.3, 4.6, 3.9, 5.6};
  //  double N[] = {22, 10, 16, 11};
  //  double SD[] = {2.11, 3.03, 2.03, 1.85};
  //  isIncreasing = false;

  // c80
  //   double D[] = {0, 25, 100, 400};
  //   double Y[] = {430.6, 431.2, 426.5, 412};
  //   double N[] = {48, 47, 49, 46};
  //   double SD[] = {28.4, 25, 30, 30.6};

  // c90
  //   double D[] = {0, 125, 250, 500, 1000, 1500};
  //   double Y[] = {352.2, 350.6, 338.8, 343.5, 330.1, 312.5};
  //   double N[] = {10, 10, 10, 10, 10, 10};
  //   double SD[] = {19.9, 11.4, 20.3, 15.2, 25, 21.6};

  // c100
  double D[] = {0, 62.5, 125, 250, 500};
  double Y[] = {24.3, 27, 31.4, 39.3, 54.2};
  double N[] = {10, 10, 10, 10, 10};
  double SD[] = {4.93, 3.16, 7.05, 13.2, 25.8};

  //    //c101b
  //  double D[] = {0, 0.156, 0.312, 0.625, 1.25, 2.5};
  //  double Y[] = {65.3, 74, 77.3, 81.3, 87.5, 92.67};
  //  double N[] = {10, 10, 10, 10, 10, 9};
  //  double SD[] = {10.18253407, 9.550078534, 16.98143104, 9.834683523, 14.60972279, 8.04};
  //  isIncreasing = true;

  // c102
  //   double D[] = {0,0.156, 0.312, 0.625, 1.25, 2.5};
  //   double Y[] = {62.6, 60.44, 57.9, 63.3, 81.9, 112.57};
  //   double N[] = {10, 9, 10, 10, 10, 7};
  //   double SD[] = {10.75174404, 6.51, 4.110960958, 4.996398703, 8.28516747, 22.54180117};

  // c103
  //   double D[] = {0, 0.156, 0.312, 0.625, 1.25, 2.5};
  //   double Y[] = {136.4, 156.1, 182.8, 184.2, 281.1, 262.4};
  //   double N[] = {9, 9, 10, 10, 10, 7};
  //   double SD[] = {18.6, 24, 36.68242086, 33.20391543, 72.41615842, 60.05855476};

  // c104
  //   double D[] = {0, 0.156, 0.312, 0.625, 1.25, 2.5};
  //   double Y[] = {35.5, 39.32, 42.61, 45.56, 54.77, 67.9};
  //   double N[] = {10, 10, 10, 10, 10, 10};
  //   double SD[] = {3.06740933, 1.67600716, 1.77087549, 2.656313235, 2.150348809, 3.763110416};

  // c105b
  //  double D[] = {0, 0.156, 0.312, 0.625, 1.25, 2.5};
  //  double Y[] = {33.52, 37.66, 40.08, 44.25, 50.84, 67.75};
  //  double N[] = {10, 10, 10, 10, 10, 10};
  //  double SD[] =
  //  {2.37170824512628, 2.81442711754986, 1.77087548969429, 2.59306768133807, 2.11872603231281, 2.84604989415154};

  // c105b truncated
  //  double D[] = {0, 0.156, 0.312, 0.625, 1.25, 2.5};
  //  double Y[] = {33.52, 37.66, 40.08, 44.25, 50.84, 67.75};
  //  double N[] = {10, 10, 10, 10, 10, 10};
  //  double SD[] = {2.372, 2.814, 1.771, 2.593, 2.119, 2.846};
  //  isIncreasing = true;

  // c106
  //  double D[] = {0, 0.125, 0.25, 0.5};
  //  double Y[] = {8766, 8831, 9215, 9906};
  //  double N[] = {8, 8, 8, 8};
  //  double SD[] = {953.179941039466, 1029.54747340761, 647.709811566878, 1100.25815152627};

  // c107
  //  double D[] = {0, 0.125, 0.25, 0.5};
  //  double Y[] = {7347, 8052, 8467, 9124};
  //  double N[] = {8, 8, 8, 8};
  //  double SD[] = {664.6803743, 933.3809512, 842.8712832, 449.7199128};

  // c108
  //  double D[] = {0, 0.125, 0.25, 0.5};
  //  double Y[] = {8390, 8342, 10114, 11633};
  //  double N[] = {8, 8, 8, 8};
  //  double SD[] = {393.1513703, 639.2245302, 987.1210665, 941.8662325};

  // c109
  //  double D[] = {0, 0.125, 0.25, 0.5};
  //  double Y[] = {4.02, 4.06, 4.35, 4.68};
  //  double N[] = {8, 8, 8, 8};
  //  double SD[] = {0.282842712, 0.282842712, 0.282842712, 0.339411255};

  // c110
  //  double D[] = {0, 0.125, 0.25, 0.5};
  //  double Y[] = {3.42, 3.77, 3.86, 4.19};
  //  double N[] = {8, 8, 8, 8};
  //  double SD[] = {0.254558441, 0.282842712, 0.254558441, 0.169705627};

  // c111
  //  double D[] = {0, 0.125, 0.25, 0.5};
  //  double Y[] = {3.85, 3.94, 4.6, 5.21};
  //  double N[] = {8, 8, 8, 8};
  //  double SD[] = {0.141421356, 0.113137085, 0.367695526, 0.282842712};

  // c113
  //  double D[] = {0, 0.04464, 0.0893, 0.179, 0.36, 0.71};
  //  double Y[] = {4.83, 5.01, 5.61, 6.14, 7.31, 8.76};
  //  double N[] = {8, 8, 8, 8, 8, 8};
  //  double SD[] = {0.22627417, 0.282842712, 0.169705627, 0.311126984, 0.282842712, 0.509116882};

  // c112
  //    double D[] = {0, 0.04464, 0.0893, 0.179, 0.36, 0.71};
  //    double Y[] = {1122, 1198, 1415, 1419, 1768, 2117};
  //    double N[] = {8, 8, 8, 8, 8, 8};
  //    double SD[] = {86.45987002, 125.8650071, 152.452222, 168.0085712, 168.8570993, 211.5663489};

  //  //c114
  //  double D[] = {0, 0.25, 0.5, 1, 2, 4, 8, 16, 32};
  //  double Y[] = {6.7, 6.4, 6.9, 7.5, 7.7, 9.7, 11.5, 12.9, 13.8};
  //  double N[] = {13, 12, 10, 13, 14, 13, 11, 13, 12};
  //  double SD[] = {0.360555128, 1.039230485, 0.948683298, 0.360555128, 0.748331477,
  //  0.360555128, 1.658312395, 1.44222051, 3.117691454};

  // c115
  //  double D[] = {0, 0.03, 0.1, 0.3, 1, 3, 6.4, 12.8};
  //  double Y[] = {6.7, 6.7, 6.9, 6.8, 7.9, 10.3, 13.8, 15.2};
  //  double N[] = {12, 11, 12, 12, 14, 12, 14, 10};
  //  double SD[] = {0.692820323, 0.663324958, 0.346410162, 0.692820323,
  //  0.374165739, 1.385640646, 1.122497216, 1.264911064};

  // c116
  //  double D[] = {0, 0.44, 3.55, 48, 92.9};
  //  double Y[] = {5.2, 5.08, 5.09, 8.29, 11.5};
  //  double N[] = {12, 10, 15, 10, 6};
  //  double SD[] = {0.13, 0.19, 0.15, 0.16, 0.24};

  // c117
  //  double D[] = {0, 0.1, 0.5, 1.1};
  //  double Y[] = {43.85, 43.51, 40.04, 35.09};
  //  double N[] = {37, 35, 43, 42};
  //  double SD[] = {2.69, 2.86, 3, 2.56};

  // c118
  //  double D[] =
  //  {0,0,0,0,0,0,0,0,0,0,3.12,3.12,3.12,3.12,3.12,3.12,3.12,3.12,3.12,3.12,6.25,6.25,6.25,6.25,6.25,6.25,6.25,6.25,6.25,6.25,12.5,12.5,12.5,12.5,12.5,12.5,12.5,12.5,12.5,12.5,25,25,25,25,25,25,25,25,25,25,50,50,50,50,50,50,50,50,50,50};
  //  double Y[] =
  //  {0.9805,1.8726,1.2946,1.4332,1.8938,1.5495,1.0806,1.758,1.5236,1.835,0.9696,1.1148,1.4757,1.7458,1.0309,1.6299,1.7097,1.4108,1.1338,1.0123,1.2254,1.9975,1.2686,1.1283,1.8501,1.0474,1.2585,1.3154,0.6003,1.0602,0.7941,1.1935,1.1676,1.0943,1.052,1.1097,0.6617,1.0424,1.063,0.9127,1.0893,0.9427,0.8838,1.0599,0.967,0.8348,1.2608,1.1349,0.7089,1.7656,1.0003,1.2963,0.9296,0.6975,0.446,0.9864,0.7209,1.1935,1.0668,1.0383};
  //  double N[1];
  //  double SD[1];

  // BMDS-165 Assaf
  //  double D[] = {1e-3, 0.02, 0.06, 0.18, 0.54, 1.62, 4.86};
  //  double Y[] = {0, 0.0428, 0.1072, 0.1968, 0.5409, 1, 1};
  //  double N[] = {0,0,0,0,0,0,0};
  //  double SD[] = {0,0,0,0,0,0,0};
  //  isIncreasing = true;
  //  double N[1];
  //  double SD[1];

  // Exact Model fit
  //  double D[] = {0,1,2,3};
  //  double Y[] = {1,2,3,4};
  //  double N[] = {10,10,10,10};
  //  double SD[] = {0.1,0.1,0.1,0.1};

  // Allen funky dataset
  //  double D[] = {0, 0.1, 0.5, 1.1};
  //  double Y[] = {43.85, 43.51, 40.04, 35.09};
  //  double N[] = {37, 35, 43, 42};
  //  double SD[] = {2.69, 2.86, 3, 2.56};

  //  double D[] = {0,25,50};
  //  double Y[] = {7.96, 9.65, 10.07};
  //  double N[] = {10,10,10};
  //  double SD[] = {3.26,3.14,3.14};

  //  double D[] = {0,50,100};
  //  double Y[] = {7.97, 9.82, 10.34};
  //  double N[] = {10,10,10};
  //  double SD[] = {2.85,2.8,2.91};

  //  double D[] = {0,50,100,200};
  //  double Y[] = {7.95, 7.6, 9.4, 9.06};
  //  double N[] = {10,10,10,10};
  //  double SD[] = {2.89,2.56,2.5,2.62};

  //  double D[] = {0, 50, 400};
  //  double Y[] = {5.26, 5.76, 9.23};
  //  double N[] = {20, 20, 20};
  //  double SD[] = {2.23, 1.47, 1.56};
  /////////////////////////////////////////////////
  // END USER INPUT
  ///////////////////////////////////////////////////

  struct continuous_analysis anal;
  int numDataRows = sizeof(D) / sizeof(D[0]);

  if (!detectAdvDir) {
    anal.isIncreasing = isIncreasing;
  }

  // check data array sizes for consistency
  size_t numElementsY = sizeof(Y) / sizeof(Y[0]);
  if (suffStat) {
    size_t numElementsN = sizeof(N) / sizeof(N[0]);
    size_t numElementsSD = sizeof(SD) / sizeof(SD[0]);
    if (numDataRows != numElementsY || numElementsY != numElementsN ||
        numElementsN != numElementsSD) {
      printf("Number of data elements are not consistent\nExiting Code\n");
      exit(-1);
    }
  } else {
    if (numDataRows != numElementsY) {
      printf("Number of data elements are not consistent\nExiting Code\n");
      exit(-1);
    }
  }

  // priors defined columnwise
  int prCols = 5;

  // define priors/parameter constraints
  int numParms;
  printf("model = %d\n", model);
  switch (model) {
    case hill:
      numParms = 6;
      break;
    case exp_3:
      // numParms = 5;
      // break;
    case exp_5:
      numParms = 6;
      break;
    case power:
      numParms = 5;
      break;
    case funl:
      numParms = 0;  // FIX THIS
      break;
    case polynomial:
      numParms = 3 + degree;
      break;
    default:
      printf("error in numParms\n");
      return;
  }
  if (dist == normal || dist == log_normal) {
    numParms -= 1;
  }

  printf("numParms = %d\n", numParms);
  double *pr;

  //  //BAYESIAN
  //  //hill;
  //  double prBayesianHill[] =
  //  {2,1,2,2,2,1,0,1,-0.69315,0.405465,0,0,1,2,1,0.2501,1,1,0,-18,0,0,0,-18,18,18,18,18,18,18};
  //  double prBayesianHillNCV[] =
  //  {2,1,2,2,1,0,1,-0.69315,0.405465,0,1,2,1,0.2501,1,0,-18,0,0,-18,18,18,18,18,18};
  //  //Power
  //  double prBayesianPower[] = {2,1,2,1,0,0,0.405465,0,1,1,0.5,1,0,-10000,0,-18,1e6,1e4,40,18};
  //  double prBayesianPowerNCV[] =
  //  {2,1,2,2,1,0,0,0.405465,0,0,1,1,0.5,0.2501,1,0,-10000,0,0,-18,1e6,1e4,40,18,18};
  //  //double pr[] = {2,1,2,1,0,0,0,0,0.1,1,0.5,2,0,-1e2,0,-18,100,1e2,40,18};  //Matt
  //  //funl
  //  //Poly
  //  double prBayesianPoly1[] = {2,1,1,0,0,0,1,2,1,0,-10000,-18,1e6,1e4,18};
  //  double prBayesianPoly2[] = {2,1,1,1,0,0,0,0,1,2,2,1,0,-10000,-10000,-18,1e6,1e4,1e4,18};
  //  //poly 2 double prBayesianPoly3[] =
  //  {2,1,1,1,1,0,0,0,0,0,1,2,2,2,1,0,-10000,-10000,-10000,-18,1e6,1e4,1e4,1e4,18}; //poly 3 double
  //  prBayesianPoly4[] =
  //  {2,1,1,1,1,1,0,0,0,0,0,0,1,2,2,2,1,1,0,-10000,-10000,-10000,-10000,-18,1e6,1e4,1e4,1e4,1e4,18};
  //  //poly 4 double prBayesianPoly5[] =
  //  {2,1,1,1,1,1,1,0,0,0,0,0,0,0,1,2,2,2,2,1,1,0,-10000,-10000,-10000,-10000,-10000,-18,1e6,1e4,1e4,1e4,1e4,1e4,18};
  //  //poly 5 double prBayesianPoly1NCV[] =
  //  {2,1,2,1,0,0,0,0,1,2,0.2501,1,0,-10000,0,-18,1e6,1e4,18,18}; double prBayesianPoly2NCV[] =
  //  {2,1,1,2,1,0,0,0,0,0,1,2,2,0.2501,1,0,-10000,-10000,0,-18,1e6,1e4,1e4,18,18}; //poly 2 double
  //  prBayesianPoly3NCV[] =
  //  {2,1,1,1,2,1,0,0,0,0,0,0,1,2,2,2,0.2501,1,0,-10000,-10000,-10000,0,-18,1e6,1e4,1e4,1e4,18,18};
  //  //poly 3 double prBayesianPoly4NCV[] =
  //  {2,1,1,1,1,2,1,0,0,0,0,0,0,0,1,2,2,2,2,0.2501,1,0,-10000,-10000,-10000,-10000,0,-18,1e6,1e4,1e4,1e4,1e4,18,18};
  //  //poly 4 double prBayesianPoly5NCV[] =
  //  {2,1,1,1,1,1,2,1,0,0,0,0,0,0,0,0,1,2,2,2,2,2,0.2501,1,0,-10000,-10000,-10000,-10000,-10000,0,-18,1e6,1e4,1e4,1e4,1e4,1e4,18,18};
  //  //poly 5
  //  //EXP3 & 5
  //  double prBayesianExp5[] =
  //  {2,2,1,2,1,0,0,0,0,0,1,1,1,0.2501,0.5,1,0,0,-20,0,-18,1e6,100,20,18,18}; double
  //  prBayesianExp5NCV[] =
  //  {2,2,1,2,2,1,0,0,0,0,0,0,1,1,1,0.2501,0.5,1,0,0,-20,0,0,-18,1e6,100,20,18,18,18};
  //  //UNRESTRICTED FREQ
  //  //Hill
  //  double prUFreqHillNormal[] =
  //  {0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,-1e8,-1000,0,1e-8,-1e3,1e8,1000,30,18,1000};  //normal dist
  //  double prUFreqHillNormalNCV[] =
  //  {0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,1,-1e8,-1e8,0,1e-8,-1e3,-1e3,1e8,1e8,30,18,1000,1000};
  //  //normal dist double prUFreqHillLognormal[] =
  //  {0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,1e-8,-1e8,0,1e-8,-1e3,1e8,1e8,100,100,1000};  //normal dist
  //  //Power
  //  double prUFreqPower[] = {0,0,0,0,0,0,0,0,1,1,1,1,1e-8,-1e8,1e-8,-1000,1e8,1e8,100,1000};
  //  double prUFreqPowerNCV[] =
  //  {0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1e-8,-1e8,1e-8,-1000,-1000,1e8,1e8,100,1000,1000};
  //  //funl;
  //  //priors for auto detect adv dir
  //  double prUFreqPoly1[] = {0,0,0,0,0,0,1,1,1,-1e8,-1e8,-1e8,1e8,1e8,1e8}; //poly 1
  //  double prUFreqPoly2[] = {0,0,0,0,0,0,0,0,1,1,1,1,-1e-8,-1e8,-1e8,-1e8,1e8,1e8,1e8,1e8}; //poly
  //  2 double prUFreqPoly3[] =
  //  {0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,-1e8,-1e8,-1e8,-1e8,-1e8,1e8,1e8,1e8,1e8,1e8}; //poly 3 double
  //  prUFreqPoly4[] =
  //  {0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,-1e8,-1e8,-1e8,-1e8,-1e8,-1e8,1e8,1e8,1e8,1e8,1e8,1e8};
  //  //poly 4 double prUFreqPoly5[] =
  //  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,-1e8,-1e8,-1e8,-1e8,-1e8,-1e8,-1e8,1e8,1e8,1e8,1e8,1e8,1e8,1e8};
  //  //poly 5 double prUFreqPoly1NCV[] =
  //  {0,0,0,0,0,0,0,0,1,1,1,1,-1e8,-1e8,-1000,-1e8,1e8,1e8,1000,1e8}; //poly 1 double
  //  prUFreqPoly2NCV[] =
  //  {0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,-1e8,-1e8,-1e8,-1000,-1e8,1e8,1e8,1e8,1000,1e8}; //poly 2
  //  double prUFreqPoly3NCV[] =
  //  {0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,-1e8,-1e8,-1e8,-1e8,-1000,-1e8,1e8,1e8,1e8,1e8,1000,1e8};
  //  //poly 3 double prUFreqPoly4NCV[] =
  //  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,-1e8,-1e8,-1e8,-1e8,-1e8,-1000,-1e8,1e8,1e8,1e8,1e8,1e8,1000,1e8};
  //  //poly 4 double prUFreqPoly5NCV[] =
  //  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,-1e8,-1e8,-1e8,-1e8,-1e8,-1e8,-1000,-1e8,1e8,1e8,1e8,1e8,1e8,1e8,1000,1e8};
  //  //poly 5
  //
  //  //RESTRICTED FREQ
  //  //EXP3
  //  //EXP5
  //  double prRFreqExp5Normal[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-18,1,-18,1e6,100,18,18,18};
  //  double
  //  prRFreqExp5NormalNCV[]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-18,1,-18,-18,1e6,100,18,18,18,18};
  //  //HILL
  //  double prRFreqHillNormal[] =
  //  {0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,-1e8,-1e8,0,1,-1e3,1e8,1e8,30,18,1000};  //normal dist double
  //  prRFreqHillNormalNCV[] =
  //  {0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,1,-1e8,-1000,0,1,-1e3,-1e3,1e8,1000,30,18,1000,1000};
  //  //normal dist double prRFreqHillLognormal[] =
  //  {0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,1e-8,-1e8,0,1,-1e3,1e8,1e8,100,100,1000};  //normal dist
  //  //POWER
  //  double prRFreqPower[] = {0,0,0,0,0,0,0,0,1,1,1,1,1e-8,-1e8,1,-1000,1e8,1e8,100,1000};  //SEG
  //  FAULT double prRFreqPowerNCV[] =
  //  {0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1e-8,-1e8,1,-1000,-1000,1e8,1e8,100,1000,1000};  //SEG FAULT
  //  //funl;
  //  //POLY
  //  //priors for auto detect adv dir
  //  double prRFreqPoly1[] = {0,0,0,0,0,0,1,1,1,-1e8,-1e8,-1e8,1e8,1e8,1e8}; //poly 1
  //  //double prRFreqPoly2[] = {0,0,0,0,0,0,0,0,1,1,1,1,-1e-8,-1e8,-1000,-1e8,1e8,1e8,1000,1e8};
  //  //poly 2 double prRFreqPoly2[] =
  //  {0,0,0,0,0,0,0,0,1,1,1,1,-1e8,-1e8,-1e8,-1e8,1e8,1e8,1e8,1e8}; //poly 2 double prRFreqPoly3[]
  //  = {0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,-1e8,-1e8,-1e8,-1e8,-1e8,1e8,1e8,1e8,1e8,1e8}; //poly 3
  //  double prRFreqPoly4[] =
  //  {0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,-1e8,-1e8,-1e8,-1e8,-1e8,-1e8,1e8,1e8,1e8,1e8,1e8,1e8};
  //  //poly 4 double prRFreqPoly5[] =
  //  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,-1e8,-1e8,-1e8,-1e8,-1e8,-1e8,-1e8,1e8,1e8,1e8,1e8,1e8,1e8,1e8};
  //  //poly 5 double prRFreqPoly1NCV[] =
  //  {0,0,0,0,0,0,0,0,1,1,1,1,-1e8,-1e8,1000,-1e8,1e8,1e8,1000,1e8}; //poly 1 double
  //  prRFreqPoly2NCV[] =
  //  {0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,-1e8,-1e8,-1e8,-1000,-1e8,1e8,1e8,1e8,1000,1e8}; //poly 2
  //  double prRFreqPoly3NCV[] =
  //  {0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,-1e8,-1e8,-1e8,-1e8,-1000,-1e8,1e8,1e8,1e8,1e8,1000,1e8};
  //  //poly 3 double prRFreqPoly4NCV[] =
  //  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,-1e8,-1e8,-1e8,-1e8,-1e8,-1000,-1e8,1e8,1e8,1e8,1e8,1e8,1000,1e8};
  //  //poly 4 double prRFreqPoly5NCV[] =
  //  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,-1e8,-1e8,-1e8,-1e8,-1e8,-1e8,-1000,-1e8,1e8,1e8,1e8,1e8,1e8,1e8,1000,1e8};
  //  //poly 4

  printf("starting priors\n");

  if (modelType == 1) {
    // frequentist
    if (restricted) {
      printf("choosing frequentist restricted priors\n");
      switch (model) {
        case hill:
          anal.model = hill;
          if (dist == normal || dist == log_normal) {
            // normal
            anal.prior = prRFreqHillNormal;
          } else {
            //} else if (dist == normal_ncv){
            // normal NCV
            anal.prior = prRFreqHillNormalNCV;
            //} else {
            //  //lognormal
            //  anal.prior = prRFreqHillLognormal;
          }
          break;
        case exp_3:
          anal.model = exp_3;
          //          if (dist == normal || dist == log_normal){
          //            anal.prior = prRFreqExp5Normal;
          //          } else {
          //            anal.prior = prRFreqExp5NormalNCV;
          //          }
          if (dist == normal) {
            anal.prior = prRFreqExp5Normal;
          } else if (dist == normal_ncv) {
            anal.prior = prRFreqExp5NormalNCV;
          } else {
            anal.prior = prRFreqExp5Lognormal;
          }
          break;
        case exp_5:
          anal.model = exp_5;
          //          if (dist == normal || dist == log_normal){
          //            anal.prior = prRFreqExp5Normal;
          //          } else {
          //            anal.prior = prRFreqExp5NormalNCV;
          //          }
          if (dist == normal) {
            anal.prior = prRFreqExp5Normal;
          } else if (dist == normal_ncv) {
            anal.prior = prRFreqExp5NormalNCV;
          } else {
            anal.prior = prRFreqExp5Lognormal;
          }
          break;
        case power:
          anal.model = power;
          if (dist == normal || dist == log_normal) {
            anal.prior = prRFreqPower;
          } else {
            anal.prior = prRFreqPowerNCV;
          }
          break;
        case funl:
          break;
        case polynomial:
          printf("choosing polynomial model\n");
          anal.model = polynomial;
          anal.degree = degree;
          if (detectAdvDir) {
            if (dist == normal || dist == log_normal) {
              printf("using advDir auto normal or log_normal dist priors\n");
              if (degree == 1) {
                anal.prior = prRFreqPoly1;
              } else if (degree == 2) {
                anal.prior = prRFreqPoly2;
              } else if (degree == 3) {
                anal.prior = prRFreqPoly3;
              } else if (degree == 4) {
                anal.prior = prRFreqPoly4;
              } else if (degree == 5) {
                anal.prior = prRFreqPoly5;
              } else {
                printf("poly restricted normal/lognormal degree error\n");
                return;
              }
            } else {
              printf("using advDir auto normal_ncv dist priors\n");
              if (degree == 1) {
                anal.prior = prRFreqPoly1NCV;
              } else if (degree == 2) {
                anal.prior = prRFreqPoly2NCV;
              } else if (degree == 3) {
                anal.prior = prRFreqPoly3NCV;
              } else if (degree == 4) {
                anal.prior = prRFreqPoly4NCV;
              } else if (degree == 5) {
                anal.prior = prRFreqPoly5NCV;
              } else {
                printf("poly restricted normal NCV degree error\n");
                return;
              }
            }
          } else {
            if (anal.isIncreasing) {
              if (dist == normal || dist == log_normal) {
                printf("using advDir up normal or log_normal dist priors\n");
                if (degree == 1) {
                  anal.prior = prRFreqPoly1Up;
                } else if (degree == 2) {
                  anal.prior = prRFreqPoly2Up;
                } else if (degree == 3) {
                  anal.prior = prRFreqPoly3Up;
                } else if (degree == 4) {
                  anal.prior = prRFreqPoly4Up;
                } else if (degree == 5) {
                  anal.prior = prRFreqPoly5Up;
                } else {
                  printf("poly restricted normal/lognormal degree error\n");
                  return;
                }
              } else {
                printf("using advDir up normal_ncv dist priors\n");
                if (degree == 1) {
                  anal.prior = prRFreqPoly1NCVUp;
                } else if (degree == 2) {
                  anal.prior = prRFreqPoly2NCVUp;
                } else if (degree == 3) {
                  anal.prior = prRFreqPoly3NCVUp;
                } else if (degree == 4) {
                  anal.prior = prRFreqPoly4NCVUp;
                } else if (degree == 5) {
                  anal.prior = prRFreqPoly5NCVUp;
                } else {
                  printf("poly restricted normal NCV degree error\n");
                  return;
                }
              }

            } else {
              if (dist == normal || dist == log_normal) {
                printf("using advDir down normal or log_normal dist priors\n");
                if (degree == 1) {
                  anal.prior = prRFreqPoly1Down;
                } else if (degree == 2) {
                  anal.prior = prRFreqPoly2Down;
                } else if (degree == 3) {
                  printf("using prRFreqPoly3Down\n");
                  anal.prior = prRFreqPoly3Down;
                } else if (degree == 4) {
                  anal.prior = prRFreqPoly4Down;
                } else if (degree == 5) {
                  anal.prior = prRFreqPoly5Down;
                } else {
                  printf("poly restricted normal/lognormal degree error\n");
                  return;
                }
              } else {
                printf("using advDir down normal_ncv dist priors\n");
                if (degree == 1) {
                  anal.prior = prRFreqPoly1NCVDown;
                } else if (degree == 2) {
                  anal.prior = prRFreqPoly2NCVDown;
                } else if (degree == 3) {
                  anal.prior = prRFreqPoly3NCVDown;
                } else if (degree == 4) {
                  anal.prior = prRFreqPoly4NCVDown;
                } else if (degree == 5) {
                  anal.prior = prRFreqPoly5NCVDown;
                } else {
                  printf("poly restricted normal NCV degree error\n");
                  return;
                }
              }
            }
          }
          break;
        default:
          printf("error with restricted models\n");
          return;
      }
    } else {
      // unrestricted
      switch (model) {
        case hill:
          anal.model = hill;
          if (dist == normal) {
            // normal
            anal.prior = prUFreqHillNormal;
          } else if (dist == normal_ncv) {
            // normal NCV
            anal.prior = prUFreqHillNormalNCV;
          } else {
            // lognormal
            anal.prior = prUFreqHillLognormal;
          }
          break;
        case exp_3:
          printf("cannot run unrestricted exponential models\n");
          return;
          // break;
        case exp_5:
          printf("cannot run unrestricted exponential models\n");
          return;
          // break;
        case power:
          anal.model = power;
          if (dist == normal || dist == log_normal) {
            anal.prior = prUFreqPower;
          } else {
            anal.prior = prUFreqPowerNCV;
          }
          break;

        case funl:
          break;
        case polynomial:
          printf("choosing polynomial model\n");
          anal.model = polynomial;
          anal.degree = degree;
          // if (detectAdvDir){
          if (dist == normal || dist == log_normal) {
            printf("prior with normal or lognormal dist\n");
            if (degree == 1) {
              anal.prior = prUFreqPoly1;
            } else if (degree == 2) {
              anal.prior = prUFreqPoly2;
            } else if (degree == 3) {
              anal.prior = prUFreqPoly3;
            } else if (degree == 4) {
              anal.prior = prUFreqPoly4;
            } else if (degree == 5) {
              anal.prior = prUFreqPoly5;
            } else {
              printf("poly unrestricted normal/lognormal degree error\n");
              return;
            }
          } else {
            if (degree == 1) {
              anal.prior = prUFreqPoly1NCV;
            } else if (degree == 2) {
              anal.prior = prUFreqPoly2NCV;
            } else if (degree == 3) {
              anal.prior = prUFreqPoly3NCV;
            } else if (degree == 4) {
              anal.prior = prUFreqPoly4NCV;
            } else if (degree == 5) {
              anal.prior = prUFreqPoly5NCV;
            } else {
              printf("poly restricted normal NCV degree error\n");
              return;
            }
          }
          //}
          break;

        default:
          printf("error with unrestricted model\n");
          return;
      }
    }
  } else {
    // bayesian
    switch (model) {
      case hill:
        anal.model = hill;
        if (dist == normal || dist == log_normal) {
          // normal
          anal.prior = prBayesianHill;
        } else {
          // normal NCV
          anal.prior = prBayesianHillNCV;
        }
        break;
      case exp_3:
        anal.model = exp_3;
        if (dist == normal || dist == log_normal) {
          // normal
          anal.prior = prBayesianExp5;
        } else {
          // normal NCV
          anal.prior = prBayesianExp5NCV;
        }
        break;
      case exp_5:
        anal.model = exp_5;
        if (dist == normal || dist == log_normal) {
          // normal
          anal.prior = prBayesianExp5;
        } else {
          // normal NCV
          anal.prior = prBayesianExp5NCV;
        }
        break;
      case power:
        anal.model = power;
        if (dist == normal || dist == log_normal) {
          // normal
          anal.prior = prBayesianPower;
        } else {
          // normal NCV
          anal.prior = prBayesianPowerNCV;
        }
        break;
      case funl:
        anal.model = funl;
        printf("FUNL model has not been implemented in BMDS\n");
        break;
      case polynomial:
        anal.model = polynomial;
        anal.degree = degree;
        if (dist == normal || dist == log_normal) {
          // normal
          printf("using Bayesian normal or log_normal dist priors\n");
          if (degree == 1) {
            anal.prior = prBayesianPoly1;
          } else if (degree == 2) {
            anal.prior = prBayesianPoly2;
          } else if (degree == 3) {
            anal.prior = prBayesianPoly3;
          } else if (degree == 4) {
            anal.prior = prBayesianPoly4;
          } else if (degree == 5) {
            anal.prior = prBayesianPoly5;
          } else {
            printf("poly restricted normal/lognormal degree error\n");
            return;
          }
        } else {
          // normal NCV
          printf("using Bayesian normal_ncv dist priors\n");
          if (degree == 1) {
            anal.prior = prBayesianPoly1NCV;
          } else if (degree == 2) {
            anal.prior = prBayesianPoly2NCV;
          } else if (degree == 3) {
            anal.prior = prBayesianPoly3NCV;
          } else if (degree == 4) {
            anal.prior = prBayesianPoly4NCV;
          } else if (degree == 5) {
            anal.prior = prBayesianPoly5NCV;
          } else {
            printf("poly restricted normal/lognormal degree error\n");
            return;
          }
        }
        break;
    }
  }

  //  printf("initial priors\n");
  //  for (int i=0; i<numParms * anal.prior_cols; i++){
  //    printf("%f,",anal.prior[i]);
  //  }

  //  printf("finished with priors\n");

  printf("prior b4 adj:\n");
  for (int i = 0; i < prCols * numParms; i++) {
    printf("%.9f\n", anal.prior[i]);
  }

  // parms array declared
  //  int numParms = sizeof(pr)/sizeof(pr[0])/prCols;
  // double parms[numParms];
  double *parms = new double[numParms];

  // declare analysis
  anal.Y = Y;
  anal.n = numDataRows;
  if (suffStat) {
    anal.n_group = N;
    anal.sd = SD;
  }
  anal.doses = D;
  anal.disttype = dist;
  if (!detectAdvDir) {
    anal.isIncreasing = isIncreasing;
  }
  anal.alpha = alpha;
  anal.BMD_type = BMD_type;  // 1=absdev, 2 = stddev, 3 = reldev, 4 = pt, 5 = extra, 6 =
                             // hybrid_extra, 7 = hybrid_added   from src/include/cmodeldefs.h
  anal.BMR = BMRF;
  anal.samples = 0;  // num MCMC samples
  anal.tail_prob = 0.01;
  anal.suff_stat = suffStat;
  anal.parms = numParms;
  anal.prior_cols = prCols;
  anal.transform_dose = 0;

  struct continuous_model_result res;
  res.model = anal.model;
  res.nparms = anal.parms;
  res.parms = parms;
  res.dist_numE = 100;

  // double cov[numParms*numParms];
  // double bmd_dist[res.dist_numE*2];
  double *cov = new double[numParms * numParms];
  double *bmd_dist = new double[res.dist_numE * 2];
  res.cov = cov;
  res.bmd_dist = bmd_dist;

  struct BMDS_results BMDSres;
  // bool bounded[anal.parms];
  // double stdErr[anal.parms];
  // double lowerConf[anal.parms];
  // double upperConf[anal.parms];
  // bool* bounded = new bool[anal.parms];
  // double* stdErr = new double[anal.parms];
  // double* lowerConf = new double[anal.parms];
  // double* upperConf = new double[anal.parms];
  // set all parms as unbounded initially
  for (int i = 0; i < anal.parms; i++) {
    // bounded[i] = false
    BMDSres.bounded.push_back(false);
    ;
    // stdErr[i] = -9999.0;
    BMDSres.stdErr.push_back(BMDS_MISSING);
    BMDSres.lowerConf.push_back(BMDS_MISSING);
    BMDSres.upperConf.push_back(BMDS_MISSING);
    // lowerConf[i] = -9999.0;
    // upperConf[i] = -9999.0;
  }
  // BMDSres.bounded = bounded;
  // BMDSres.stdErr = stdErr;
  // BMDSres.lowerConf = lowerConf;
  // BMDSres.upperConf = upperConf;
  BMDSres.BMD = -9999.0;
  BMDSres.BMDU = -9999.0;
  BMDSres.BMDL = -9999.0;
  BMDSres.AIC = -9999.0;

  struct continuous_GOF gof;
  int nGOF;
  if (anal.suff_stat) {
    nGOF = anal.n;
  } else {
    // determine number of unique dose groups
    nGOF = 1;
    for (int i = 1; i < anal.n; i++) {
      int j = 0;
      for (j = 0; j < i; j++) {
        if (anal.doses[i] == anal.doses[j]) break;
      }
      if (i == j) nGOF++;
    }
  }

  printf("nGOF = %d\n", nGOF);

  // double* doseGOF = new double[nGOF];
  // double* sizeGOF = new double[nGOF];
  // double* estMeanGOF = new double[nGOF];
  // double* calcMeanGOF = new double[nGOF];
  // double* obsMeanGOF = new double[nGOF];
  // double* estSDGOF = new double[nGOF];
  // double* calcSDGOF = new double[nGOF];
  // double* obsSDGOF = new double[nGOF];
  // double* resGOF = new double[nGOF];
  // double* ebLower = new double[nGOF];
  // double* ebUpper = new double[nGOF];

  // gof.dose = doseGOF;
  // gof.size = sizeGOF;
  // gof.estMean = estMeanGOF;
  // gof.calcMean = calcMeanGOF;
  // gof.obsMean = obsMeanGOF;
  // gof.estSD = estSDGOF;
  // gof.calcSD = calcSDGOF;
  // gof.obsSD = obsSDGOF;
  // gof.res = resGOF;
  // gof.n = nGOF;
  // gof.ebLower = ebLower;
  // gof.ebUpper = ebUpper;

  struct continuous_AOD aod;
  // double LL[5];
  // int nParms[5];
  // double AIC[5];
  std::vector<double> LL(5);
  std::vector<int> nParms(5);
  std::vector<double> AIC(5);
  double addConst;  // = 22.2;
  aod.LL = LL;
  aod.nParms = nParms;
  aod.AIC = AIC;
  aod.addConst = addConst;

  // double llRatio[4];
  // double DF[4];
  // double pVal[4];

  std::vector<double> llRatio(4);
  std::vector<double> DF(4);
  std::vector<double> pVal(4);
  struct testsOfInterest TOI;

  TOI.llRatio = llRatio;
  TOI.DF = DF;
  TOI.pVal = pVal;
  aod.TOI = TOI;

  printf("\n\n-----------INPUT---------------\n");
  printf("priors sent to model code:\n");
  for (int i = 0; i < prCols * numParms; i++) {
    printf("%.20f\n", anal.prior[i]);
  }
  printf("\n\ncontinuous_analysis values\n");
  printf("CA.n = %d\n", anal.n);
  printf("CA.suff_stat = %s\n", (anal.suff_stat ? "true" : "false"));
  //  printf("CA data arrays (dose, Y, n_group, sd)\n");
  //  for (int i=0; i<anal.n;i++){
  //    printf("%f, %f, %f, %f\n", anal.doses[i], anal.Y[i], anal.n_group[i], anal.sd[i]);
  //  }
  printf("CA.BMD_type = %d\n", anal.BMD_type);
  printf("CA.isIncreasing = %s\n", (anal.isIncreasing ? "true" : "false"));
  printf("CA.BMR = %.20f\n", anal.BMR);
  printf("CA.disttype = %d\n", anal.disttype);
  printf("CA.alpha = %.20f\n", anal.alpha);
  printf("CA.degree = %d\n", anal.degree);
  printf("CA.parms = %d\n", anal.parms);
  printf("CA.prior_cols = %d\n", anal.prior_cols);

  printf("\n\nData\n");
  printf("Dose, N, Mean, Std. Dev.\n");
  for (int i = 0; i < anal.n; i++) {
    printf("%.20f, %.20f, %.20f, %.20f\n", anal.doses[i], anal.n_group[i], anal.Y[i], anal.sd[i]);
  }

  printf("\n\n");
  printf("calling runBMDSContAnalysis\n");
  runBMDSContAnalysis(
      &anal, &res, &BMDSres, &aod, &gof, &detectAdvDir, &restricted, &countAllParmsOnBoundary
  );

  if (detectAdvDir) {
    printf("auto adverse direction: %s\n", anal.isIncreasing ? "increasing" : "decreasing");
  }

  printf("\n\n");
  printf("prior after adj by model code:\n");
  for (int i = 0; i < prCols * numParms; i++) {
    printf("%.20f\n", anal.prior[i]);
  }

  printf("\n\n----------OUTPUT-----------\n");
  printf("tlink BMDSres.validResult = %s\n", BMDSres.validResult ? "valid" : "invalid");
  if (BMDSres.validResult || showResultsOverride) {
    printf("\nBenchmark Dose\n");
    printf("max:  %f\n", res.max);
    printf("BMD:  %f\n", BMDSres.BMD);
    printf("Matt's BMD:  %f\n", res.bmd);
    printf("BMDL: %f\n", BMDSres.BMDL);
    printf("BMDU: %f\n", BMDSres.BMDU);
    printf("AIC:  %f\n", BMDSres.AIC);
    printf("LPP: %f\n", BMDSres.BIC_equiv);
    printf("Test 4 P-value: %f\n", aod.TOI.pVal[3]);
    printf("DOF: %f\n", aod.TOI.DF[3]);
    printf("ChiSq: %f\n", BMDSres.chisq);

    printf("\nModel Parameters\n");
    printf("# of parms: %d\n", anal.parms);
    printf("parm, estimate, bounded, std.err., lower conf, upper conf\n");
    for (int i = 0; i < anal.parms; i++) {
      printf(
          "%d, %.20f, %s, %f, %f, %f\n", i, res.parms[i], BMDSres.bounded[i] ? "true" : "false",
          BMDSres.stdErr[i], BMDSres.lowerConf[i], BMDSres.upperConf[i]
      );
      //     printf("bounded %d = %s\n", i, BMDSres.bounded[i] ? "true" : "false");
    }

    printf("\nGoodness of Fit\n");
    printf("gof.n = %d\n", gof.n);
    printf("Dose, Size, EstMed, CalcMed, ObsMean, EstSD, CalcSD, ObsSD, SR\n");
    for (int i = 0; i < gof.n; i++) {
      printf(
          "%f, %f, %f, %f, %f, %f, %f, %f, %f\n", gof.dose[i], gof.size[i], gof.estMean[i],
          gof.calcMean[i], gof.obsMean[i], gof.estSD[i], gof.calcSD[i], gof.obsSD[i], gof.res[i]
      );
    }
    printf("\nError Bars\n");
    for (int i = 0; i < gof.n; i++) {
      printf("%f, %f\n", gof.ebLower[i], gof.ebUpper[i]);
    }

    printf("\nLikelihoods of Interest\n");
    for (int i = 0; i < 5; i++) {
      printf("i:%d, LL:%f, nParms:%d, AIC:%f\n", i, aod.LL[i], aod.nParms[i], aod.AIC[i]);
    }
    printf("additive constant:%f\n", aod.addConst);

    printf("\nTests of Interest:\n");
    for (int i = 0; i < 4; i++) {
      printf(
          "i:%d, llRatio:%f, DF:%f, pVal:%f\n", i, aod.TOI.llRatio[i], aod.TOI.DF[i],
          aod.TOI.pVal[i]
      );
    }

    printf("\nBMD Dist:\n");
    for (int i = 0; i < res.dist_numE; i++) {
      printf("i:%d, perc:%f, dist:%f\n", i, res.bmd_dist[i + res.dist_numE], res.bmd_dist[i]);
    }

  } else {
    printf("Model was not run\n");
  }

  // debugContAnal(&anal);
}

void runPythonContAnalysis() {
  printf("Running python continuous analysis\n");

  // char * bmdsVersion[32];

  // version(bmdsVersion);
  string bmdsVersion = version();
  // printf("Version: %s\n", bmdsVersion);
  std::cout << "Version: " << bmdsVersion << std::endl;

  bool isIncreasing;

  ///////////////////////////////
  // USER INPUT
  //////////////////////////////

  enum cont_model model = polynomial;  // hill, exp_3, exp_5, power, funl, polynomial
  int modelType = 1;                   // 1 = frequentist, 2 = bayesian
  bool restricted = true;              // only used for frequentist models
  enum distribution dist = normal;     // normal, normal_ncv, log_normal
  bool detectAdvDir = true;            // if false then need to set isIncreasing
  bool countAllParmsOnBoundary = true;
  // isIncreasing = true;

  int degree = 2;  // for polynomial only

  double alpha = 0.05;
  double BMRF = 1.0;  // 1.0;
  int BMD_type = 2;   // 1=absdev, 2 = stddev, 3 = reldev, 4 = pt, 5 = extra, 6 = hybrid_extra, 7 =
                      // hybrid_added   from src/include/cmodeldefs.h
  ////////////////////////////////////////////////
  // cont data - suff stat: dose, Y, N, SD
  // cont data - individual: dose, response
  /////////////////////////////////////////////
  bool suffStat = true;

  // continuous1.dax
  //  double D[] = {0,25,50, 100, 200};
  //  double Y[] = {6.0, 5.2, 2.4, 1.1, 0.75};
  //  double N[] = {20, 20, 19, 20, 20};
  //  double SD[] = {1.2, 1.1, 0.81, 0.74, 0.66};
  // isIncreasing = false;

  // continuous2.dax
  // double D[] = {0,  0,  0,  0,  18, 18, 18, 18, 18, 20, 20, 20, 20,
  //               30, 30, 30, 30, 35, 35, 35, 35, 40, 40, 40, 40, 40};
  // double Y[] = {39,   38.4, 36.3, 37.1, 40.2, 45.3, 42.1, 38.3, 35.9, 42.5, 45.2, 40.1, 39.8,
  //               50.1, 53.4, 48.2, 52.1, 56.1, 50.4, 53.2, 55.2, 55.1, 59.1, 56.3, 52.9, 53.7};
  // double N[1];
  // double SD[1];
  // isIncreasing = true;

  // continuous3.dax
  //  double D[] = {0,35,105,316,625};
  //  double Y[] = {1.61,1.66,1.75,1.81,1.89};
  //  double N[] = {10,10,10,10,10};
  //  double SD[] = {0.12,0.13,0.11,0.15,0.13};

  // other test datasets
  //  double D[] = {0,50, 100, 150, 200};
  //  double Y[] = {10, 20 , 30, 40 ,50};
  //  double N[] = {100, 100, 100, 100, 100};
  //  double SD[] = {3, 4, 5, 6, 7};
  //  isIncreasing = true;

  //  double D[] = {0,50, 100, 150, 200};
  //  double Y[] = {10, 0 , -10, -20 ,-30};
  //  double N[] = {100, 100, 100, 100, 100};
  //  double SD[] = {3, 4, 5, 6, 7};
  //  isIncreasing = false;

  //  double D[] = {0,50, 100, 150, 200};
  //  double Y[] = {10, 18, 32, 38, 70};
  //  double N[] = {100, 100, 100, 100, 100};
  //  double SD[] = {3.2, 4.8, 6.5, 7.2, 8.4};
  //  isIncreasing = true;

  //  double D[] = {0,50, 100, 150, 200};
  //  double Y[] = {1, -5 , -10, -20 ,-30};
  //  double N[] = {100, 100, 100, 100, 100};
  //  double SD[] = {3, 4, 5, 6, 7};
  //  isIncreasing = false;

  // c1b
  //  double D[] = {0, 75, 250};
  //  double Y[] = {15.81, 17.91, 21.48};
  //  double N[] = {11, 11, 11};
  //  double SD[] = {2.793, 2.902, 5.771};
  //  isIncreasing = true;

  // c2
  //    double D[] = {0,75,250};
  //    double Y[] = {94.18, 91.07, 116.61};
  //    double N[] = {11,11,11};
  //    double SD[] = {14.53, 9.22, 18.31};

  // c10
  //  double D[] = {0, 75, 250};
  //  double Y[] = {3.8, 8.9, 8.9};
  //  double N[] = {22, 22, 22};
  //  double SD[] = {3.34, 5.16, 2.49};

  // c20
  //  double D[] = {0, 1, 3, 9};
  //  double Y[] = {1.037, 1.05, 1.052, 1.066};
  //  double N[] = {10, 10, 10, 10};
  //  double SD[] = {0.015, 0.01, 0.01, 0.01};

  // c40
  //  double D[] = {0, 25, 100, 400};
  //  double Y[] = {0.67, 0.68, 0.71, 0.62};
  //  double N[] = {14, 15, 15, 15};
  //  double SD[] = {0.13, 0.09, 0.11, 0.08};

  // c60
  //  double D[] = {0, 10, 50, 100, 250};
  //  double Y[] = {0.116, 0.113, 0.108, 0.108, 0.106};
  //  double N[] = {30, 30, 30, 30, 30};
  //  double SD[] = {0.006, 0.006, 0.004, 0.009, 0.008};

  // c70b
  //  double D[] = {0, 46.4, 68.1, 200};
  //  double Y[] = {6.3, 4.6, 3.9, 5.6};
  //  double N[] = {22, 10, 16, 11};
  //  double SD[] = {2.11, 3.03, 2.03, 1.85};
  //  isIncreasing = false;

  // c80
  //   double D[] = {0, 25, 100, 400};
  //   double Y[] = {430.6, 431.2, 426.5, 412};
  //   double N[] = {48, 47, 49, 46};
  //   double SD[] = {28.4, 25, 30, 30.6};

  // c90
  //   double D[] = {0, 125, 250, 500, 1000, 1500};
  //   double Y[] = {352.2, 350.6, 338.8, 343.5, 330.1, 312.5};
  //   double N[] = {10, 10, 10, 10, 10, 10};
  //   double SD[] = {19.9, 11.4, 20.3, 15.2, 25, 21.6};

  // c100
  //   double D[] = {0, 62.5, 125, 250, 500};
  //   double Y[] = {24.3, 27, 31.4, 39.3, 54.2};
  //   double N[] = {10, 10, 10, 10, 10};
  //   double SD[] = {4.93, 3.16, 7.05, 13.2, 25.8};

  //    //c101b
  //  double D[] = {0, 0.156, 0.312, 0.625, 1.25, 2.5};
  //  double Y[] = {65.3, 74, 77.3, 81.3, 87.5, 92.67};
  //  double N[] = {10, 10, 10, 10, 10, 9};
  //  double SD[] = {10.18253407, 9.550078534, 16.98143104, 9.834683523, 14.60972279, 8.04};
  //  isIncreasing = true;

  // c102
  //   double D[] = {0,0.156, 0.312, 0.625, 1.25, 2.5};
  //   double Y[] = {62.6, 60.44, 57.9, 63.3, 81.9, 112.57};
  //   double N[] = {10, 9, 10, 10, 10, 7};
  //   double SD[] = {10.75174404, 6.51, 4.110960958, 4.996398703, 8.28516747, 22.54180117};

  // c103
  //   double D[] = {0, 0.156, 0.312, 0.625, 1.25, 2.5};
  //   double Y[] = {136.4, 156.1, 182.8, 184.2, 281.1, 262.4};
  //   double N[] = {9, 9, 10, 10, 10, 7};
  //   double SD[] = {18.6, 24, 36.68242086, 33.20391543, 72.41615842, 60.05855476};

  // c104
  //   double D[] = {0, 0.156, 0.312, 0.625, 1.25, 2.5};
  //   double Y[] = {35.5, 39.32, 42.61, 45.56, 54.77, 67.9};
  //   double N[] = {10, 10, 10, 10, 10, 10};
  //   double SD[] = {3.06740933, 1.67600716, 1.77087549, 2.656313235, 2.150348809, 3.763110416};

  // c105b
  //  double D[] = {0, 0.156, 0.312, 0.625, 1.25, 2.5};
  //  double Y[] = {33.52, 37.66, 40.08, 44.25, 50.84, 67.75};
  //  double N[] = {10, 10, 10, 10, 10, 10};
  //  double SD[] =
  //  {2.37170824512628, 2.81442711754986, 1.77087548969429, 2.59306768133807, 2.11872603231281, 2.84604989415154};

  // c105b truncated
  //  double D[] = {0, 0.156, 0.312, 0.625, 1.25, 2.5};
  //  double Y[] = {33.52, 37.66, 40.08, 44.25, 50.84, 67.75};
  //  double N[] = {10, 10, 10, 10, 10, 10};
  //  double SD[] = {2.372, 2.814, 1.771, 2.593, 2.119, 2.846};
  //  isIncreasing = true;

  // c106
  //  double D[] = {0, 0.125, 0.25, 0.5};
  //  double Y[] = {8766, 8831, 9215, 9906};
  //  double N[] = {8, 8, 8, 8};
  //  double SD[] = {953.179941039466, 1029.54747340761, 647.709811566878, 1100.25815152627};

  // c107
  //  double D[] = {0, 0.125, 0.25, 0.5};
  //  double Y[] = {7347, 8052, 8467, 9124};
  //  double N[] = {8, 8, 8, 8};
  //  double SD[] = {664.6803743, 933.3809512, 842.8712832, 449.7199128};

  // c108
  //  double D[] = {0, 0.125, 0.25, 0.5};
  //  double Y[] = {8390, 8342, 10114, 11633};
  //  double N[] = {8, 8, 8, 8};
  //  double SD[] = {393.1513703, 639.2245302, 987.1210665, 941.8662325};

  // c109
  //  double D[] = {0, 0.125, 0.25, 0.5};
  //  double Y[] = {4.02, 4.06, 4.35, 4.68};
  //  double N[] = {8, 8, 8, 8};
  //  double SD[] = {0.282842712, 0.282842712, 0.282842712, 0.339411255};

  // c110
  //  double D[] = {0, 0.125, 0.25, 0.5};
  //  double Y[] = {3.42, 3.77, 3.86, 4.19};
  //  double N[] = {8, 8, 8, 8};
  //  double SD[] = {0.254558441, 0.282842712, 0.254558441, 0.169705627};

  // c111
  //  double D[] = {0, 0.125, 0.25, 0.5};
  //  double Y[] = {3.85, 3.94, 4.6, 5.21};
  //  double N[] = {8, 8, 8, 8};
  //  double SD[] = {0.141421356, 0.113137085, 0.367695526, 0.282842712};

  // c113
  //  double D[] = {0, 0.04464, 0.0893, 0.179, 0.36, 0.71};
  //  double Y[] = {4.83, 5.01, 5.61, 6.14, 7.31, 8.76};
  //  double N[] = {8, 8, 8, 8, 8, 8};
  //  double SD[] = {0.22627417, 0.282842712, 0.169705627, 0.311126984, 0.282842712, 0.509116882};

  // c112
  //    double D[] = {0, 0.04464, 0.0893, 0.179, 0.36, 0.71};
  //    double Y[] = {1122, 1198, 1415, 1419, 1768, 2117};
  //    double N[] = {8, 8, 8, 8, 8, 8};
  //    double SD[] = {86.45987002, 125.8650071, 152.452222, 168.0085712, 168.8570993, 211.5663489};

  //  //c114
  //  double D[] = {0, 0.25, 0.5, 1, 2, 4, 8, 16, 32};
  //  double Y[] = {6.7, 6.4, 6.9, 7.5, 7.7, 9.7, 11.5, 12.9, 13.8};
  //  double N[] = {13, 12, 10, 13, 14, 13, 11, 13, 12};
  //  double SD[] = {0.360555128, 1.039230485, 0.948683298, 0.360555128, 0.748331477,
  //  0.360555128, 1.658312395, 1.44222051, 3.117691454};

  // c115
  //  double D[] = {0, 0.03, 0.1, 0.3, 1, 3, 6.4, 12.8};
  //  double Y[] = {6.7, 6.7, 6.9, 6.8, 7.9, 10.3, 13.8, 15.2};
  //  double N[] = {12, 11, 12, 12, 14, 12, 14, 10};
  //  double SD[] = {0.692820323, 0.663324958, 0.346410162, 0.692820323,
  //  0.374165739, 1.385640646, 1.122497216, 1.264911064};

  // c116
  //  double D[] = {0, 0.44, 3.55, 48, 92.9};
  //  double Y[] = {5.2, 5.08, 5.09, 8.29, 11.5};
  //  double N[] = {12, 10, 15, 10, 6};
  //  double SD[] = {0.13, 0.19, 0.15, 0.16, 0.24};

  // c117
  //  double D[] = {0, 0.1, 0.5, 1.1};
  //  double Y[] = {43.85, 43.51, 40.04, 35.09};
  //  double N[] = {37, 35, 43, 42};
  //  double SD[] = {2.69, 2.86, 3, 2.56};

  // c118
  //  double D[] =
  //  {0,0,0,0,0,0,0,0,0,0,3.12,3.12,3.12,3.12,3.12,3.12,3.12,3.12,3.12,3.12,6.25,6.25,6.25,6.25,6.25,6.25,6.25,6.25,6.25,6.25,12.5,12.5,12.5,12.5,12.5,12.5,12.5,12.5,12.5,12.5,25,25,25,25,25,25,25,25,25,25,50,50,50,50,50,50,50,50,50,50};
  //  double Y[] =
  //  {0.9805,1.8726,1.2946,1.4332,1.8938,1.5495,1.0806,1.758,1.5236,1.835,0.9696,1.1148,1.4757,1.7458,1.0309,1.6299,1.7097,1.4108,1.1338,1.0123,1.2254,1.9975,1.2686,1.1283,1.8501,1.0474,1.2585,1.3154,0.6003,1.0602,0.7941,1.1935,1.1676,1.0943,1.052,1.1097,0.6617,1.0424,1.063,0.9127,1.0893,0.9427,0.8838,1.0599,0.967,0.8348,1.2608,1.1349,0.7089,1.7656,1.0003,1.2963,0.9296,0.6975,0.446,0.9864,0.7209,1.1935,1.0668,1.0383};
  //  double N[1];
  //  double SD[1];

  // BMDS-165 Assaf
  //  double D[] = {1e-3, 0.02, 0.06, 0.18, 0.54, 1.62, 4.86};
  //  double Y[] = {0, 0.0428, 0.1072, 0.1968, 0.5409, 1, 1};
  //  double N[] = {0,0,0,0,0,0,0};
  //  double SD[] = {0,0,0,0,0,0,0};
  //  isIncreasing = true;
  //  double N[1];
  //  double SD[1];

  // Exact Model fit
  //  double D[] = {0,1,2,3};
  //  double Y[] = {1,2,3,4};
  //  double N[] = {10,10,10,10};
  //  double SD[] = {0.1,0.1,0.1,0.1};

  // Allen funky dataset
  //  double D[] = {0, 0.1, 0.5, 1.1};
  //  double Y[] = {43.85, 43.51, 40.04, 35.09};
  //  double N[] = {37, 35, 43, 42};
  //  double SD[] = {2.69, 2.86, 3, 2.56};

  //  double D[] = {0,25,50};
  //  double Y[] = {7.96, 9.65, 10.07};
  //  double N[] = {10,10,10};
  //  double SD[] = {3.26,3.14,3.14};

  //  double D[] = {0,50,100};
  //  double Y[] = {7.97, 9.82, 10.34};
  //  double N[] = {10,10,10};
  //  double SD[] = {2.85,2.8,2.91};

  //  double D[] = {0,50,100,200};
  //  double Y[] = {7.95, 7.6, 9.4, 9.06};
  //  double N[] = {10,10,10,10};
  //  double SD[] = {2.89,2.56,2.5,2.62};

  //  double D[] = {0, 50, 400};
  //  double Y[] = {5.26, 5.76, 9.23};
  //  double N[] = {20, 20, 20};
  //  double SD[] = {2.23, 1.47, 1.56};

  double D[] = {0, 50, 100, 150};
  double Y[] = {10, 15, 20, 25};
  double N[] = {10, 10, 10, 10};
  double SD[] = {2, 3, 4, 5};

  /////////////////////////////////////////////////
  // END USER INPUT
  ///////////////////////////////////////////////////

  // struct continuous_analysis anal;
  struct python_continuous_analysis anal;
  anal.countAllParmsOnBoundary = countAllParmsOnBoundary;
  int numDataRows = sizeof(D) / sizeof(D[0]);

  if (!detectAdvDir) {
    anal.isIncreasing = isIncreasing;
  }

  // check data array sizes for consistency
  size_t numElementsY = sizeof(Y) / sizeof(Y[0]);
  if (suffStat) {
    size_t numElementsN = sizeof(N) / sizeof(N[0]);
    size_t numElementsSD = sizeof(SD) / sizeof(SD[0]);
    if (numDataRows != numElementsY || numElementsY != numElementsN ||
        numElementsN != numElementsSD) {
      printf("Number of data elements are not consistent\nExiting Code\n");
      exit(-1);
    }
  } else {
    if (numDataRows != numElementsY) {
      printf("Number of data elements are not consistent\nExiting Code\n");
      exit(-1);
    }
  }

  // priors defined columnwise
  int prCols = 5;

  // define priors/parameter constraints
  int numParms;
  printf("model = %d\n", model);
  switch (model) {
    case hill:
      numParms = 6;
      break;
    case exp_3:
      // numParms = 5;
      // break;
    case exp_5:
      numParms = 6;
      break;
    case power:
      numParms = 5;
      break;
    case funl:
      numParms = 0;  // FIX THIS
      break;
    case polynomial:
      numParms = 3 + degree;
      break;
    default:
      printf("error in numParms\n");
      return;
  }
  if (dist == normal || dist == log_normal) {
    numParms -= 1;
  }

  printf("numParms = %d\n", numParms);
  // double* pr;
  double *prior;

  printf("starting priors\n");

  if (modelType == 1) {
    // frequentist
    if (restricted) {
      printf("choosing frequentist restricted priors\n");
      switch (model) {
        case hill:
          anal.model = hill;
          if (dist == normal || dist == log_normal) {
            // normal
            prior = prRFreqHillNormal;
          } else {
            //} else if (dist == normal_ncv){
            // normal NCV
            prior = prRFreqHillNormalNCV;
            //} else {
            //  //lognormal
            //  anal.prior = prRFreqHillLognormal;
          }
          break;
        case exp_3:
          anal.model = exp_3;
          //          if (dist == normal || dist == log_normal){
          //            anal.prior = prRFreqExp5Normal;
          //          } else {
          //            anal.prior = prRFreqExp5NormalNCV;
          //          }
          if (dist == normal) {
            prior = prRFreqExp5Normal;
          } else if (dist == normal_ncv) {
            prior = prRFreqExp5NormalNCV;
          } else {
            prior = prRFreqExp5Lognormal;
          }
          break;
        case exp_5:
          anal.model = exp_5;
          //          if (dist == normal || dist == log_normal){
          //            anal.prior = prRFreqExp5Normal;
          //          } else {
          //            anal.prior = prRFreqExp5NormalNCV;
          //          }
          if (dist == normal) {
            prior = prRFreqExp5Normal;
          } else if (dist == normal_ncv) {
            prior = prRFreqExp5NormalNCV;
          } else {
            prior = prRFreqExp5Lognormal;
          }
          break;
        case power:
          anal.model = power;
          if (dist == normal || dist == log_normal) {
            prior = prRFreqPower;
          } else {
            prior = prRFreqPowerNCV;
          }
          break;
        case funl:
          break;
        case polynomial:
          printf("choosing polynomial model\n");
          anal.model = polynomial;
          anal.degree = degree;
          if (detectAdvDir) {
            if (dist == normal || dist == log_normal) {
              printf("using advDir auto normal or log_normal dist priors\n");
              if (degree == 1) {
                prior = prRFreqPoly1;
              } else if (degree == 2) {
                prior = prRFreqPoly2;
              } else if (degree == 3) {
                prior = prRFreqPoly3;
              } else if (degree == 4) {
                prior = prRFreqPoly4;
              } else if (degree == 5) {
                prior = prRFreqPoly5;
              } else {
                printf("poly restricted normal/lognormal degree error\n");
                return;
              }
            } else {
              printf("using advDir auto normal_ncv dist priors\n");
              if (degree == 1) {
                prior = prRFreqPoly1NCV;
              } else if (degree == 2) {
                prior = prRFreqPoly2NCV;
              } else if (degree == 3) {
                prior = prRFreqPoly3NCV;
              } else if (degree == 4) {
                prior = prRFreqPoly4NCV;
              } else if (degree == 5) {
                prior = prRFreqPoly5NCV;
              } else {
                printf("poly restricted normal NCV degree error\n");
                return;
              }
            }
          } else {
            if (anal.isIncreasing) {
              if (dist == normal || dist == log_normal) {
                printf("using advDir up normal or log_normal dist priors\n");
                if (degree == 1) {
                  prior = prRFreqPoly1Up;
                } else if (degree == 2) {
                  prior = prRFreqPoly2Up;
                } else if (degree == 3) {
                  prior = prRFreqPoly3Up;
                } else if (degree == 4) {
                  prior = prRFreqPoly4Up;
                } else if (degree == 5) {
                  prior = prRFreqPoly5Up;
                } else {
                  printf("poly restricted normal/lognormal degree error\n");
                  return;
                }
              } else {
                printf("using advDir up normal_ncv dist priors\n");
                if (degree == 1) {
                  prior = prRFreqPoly1NCVUp;
                } else if (degree == 2) {
                  prior = prRFreqPoly2NCVUp;
                } else if (degree == 3) {
                  prior = prRFreqPoly3NCVUp;
                } else if (degree == 4) {
                  prior = prRFreqPoly4NCVUp;
                } else if (degree == 5) {
                  prior = prRFreqPoly5NCVUp;
                } else {
                  printf("poly restricted normal NCV degree error\n");
                  return;
                }
              }

            } else {
              if (dist == normal || dist == log_normal) {
                printf("using advDir down normal or log_normal dist priors\n");
                if (degree == 1) {
                  prior = prRFreqPoly1Down;
                } else if (degree == 2) {
                  prior = prRFreqPoly2Down;
                } else if (degree == 3) {
                  printf("using prRFreqPoly3Down\n");
                  prior = prRFreqPoly3Down;
                } else if (degree == 4) {
                  prior = prRFreqPoly4Down;
                } else if (degree == 5) {
                  prior = prRFreqPoly5Down;
                } else {
                  printf("poly restricted normal/lognormal degree error\n");
                  return;
                }
              } else {
                printf("using advDir down normal_ncv dist priors\n");
                if (degree == 1) {
                  prior = prRFreqPoly1NCVDown;
                } else if (degree == 2) {
                  prior = prRFreqPoly2NCVDown;
                } else if (degree == 3) {
                  prior = prRFreqPoly3NCVDown;
                } else if (degree == 4) {
                  prior = prRFreqPoly4NCVDown;
                } else if (degree == 5) {
                  prior = prRFreqPoly5NCVDown;
                } else {
                  printf("poly restricted normal NCV degree error\n");
                  return;
                }
              }
            }
          }
          break;
        default:
          printf("error with restricted models\n");
          return;
      }
    } else {
      // unrestricted
      switch (model) {
        case hill:
          anal.model = hill;
          if (dist == normal) {
            // normal
            prior = prUFreqHillNormal;
          } else if (dist == normal_ncv) {
            // normal NCV
            prior = prUFreqHillNormalNCV;
          } else {
            // lognormal
            prior = prUFreqHillLognormal;
          }
          break;
        case exp_3:
          printf("cannot run unrestricted exponential models\n");
          return;
          // break;
        case exp_5:
          printf("cannot run unrestricted exponential models\n");
          return;
          // break;
        case power:
          anal.model = power;
          if (dist == normal || dist == log_normal) {
            prior = prUFreqPower;
          } else {
            prior = prUFreqPowerNCV;
          }
          break;

        case funl:
          break;
        case polynomial:
          printf("choosing polynomial model\n");
          anal.model = polynomial;
          anal.degree = degree;
          // if (detectAdvDir){
          if (dist == normal || dist == log_normal) {
            printf("prior with normal or lognormal dist\n");
            if (degree == 1) {
              prior = prUFreqPoly1;
            } else if (degree == 2) {
              prior = prUFreqPoly2;
            } else if (degree == 3) {
              prior = prUFreqPoly3;
            } else if (degree == 4) {
              prior = prUFreqPoly4;
            } else if (degree == 5) {
              prior = prUFreqPoly5;
            } else {
              printf("poly unrestricted normal/lognormal degree error\n");
              return;
            }
          } else {
            if (degree == 1) {
              prior = prUFreqPoly1NCV;
            } else if (degree == 2) {
              prior = prUFreqPoly2NCV;
            } else if (degree == 3) {
              prior = prUFreqPoly3NCV;
            } else if (degree == 4) {
              prior = prUFreqPoly4NCV;
            } else if (degree == 5) {
              prior = prUFreqPoly5NCV;
            } else {
              printf("poly restricted normal NCV degree error\n");
              return;
            }
          }
          //}
          break;

        default:
          printf("error with unrestricted model\n");
          return;
      }
    }
  } else {
    // bayesian
    switch (model) {
      case hill:
        anal.model = hill;
        if (dist == normal || dist == log_normal) {
          // normal
          prior = prBayesianHill;
        } else {
          // normal NCV
          prior = prBayesianHillNCV;
        }
        break;
      case exp_3:
        anal.model = exp_3;
        if (dist == normal || dist == log_normal) {
          // normal
          prior = prBayesianExp5;
        } else {
          // normal NCV
          prior = prBayesianExp5NCV;
        }
        break;
      case exp_5:
        anal.model = exp_5;
        if (dist == normal || dist == log_normal) {
          // normal
          prior = prBayesianExp5;
        } else {
          // normal NCV
          prior = prBayesianExp5NCV;
        }
        break;
      case power:
        anal.model = power;
        if (dist == normal || dist == log_normal) {
          // normal
          prior = prBayesianPower;
        } else {
          // normal NCV
          prior = prBayesianPowerNCV;
        }
        break;
      case funl:
        anal.model = funl;
        printf("FUNL model has not been implemented in BMDS\n");
        break;
      case polynomial:
        anal.model = polynomial;
        anal.degree = degree;
        if (dist == normal || dist == log_normal) {
          // normal
          printf("using Bayesian normal or log_normal dist priors\n");
          if (degree == 1) {
            prior = prBayesianPoly1;
          } else if (degree == 2) {
            prior = prBayesianPoly2;
          } else if (degree == 3) {
            prior = prBayesianPoly3;
          } else if (degree == 4) {
            prior = prBayesianPoly4;
          } else if (degree == 5) {
            prior = prBayesianPoly5;
          } else {
            printf("poly restricted normal/lognormal degree error\n");
            return;
          }
        } else {
          // normal NCV
          printf("using Bayesian normal_ncv dist priors\n");
          if (degree == 1) {
            prior = prBayesianPoly1NCV;
          } else if (degree == 2) {
            prior = prBayesianPoly2NCV;
          } else if (degree == 3) {
            prior = prBayesianPoly3NCV;
          } else if (degree == 4) {
            prior = prBayesianPoly4NCV;
          } else if (degree == 5) {
            prior = prBayesianPoly5NCV;
          } else {
            printf("poly restricted normal/lognormal degree error\n");
            return;
          }
        }
        break;
    }
  }

  //  printf("initial priors\n");
  //  for (int i=0; i<numParms * anal.prior_cols; i++){
  //    printf("%f,",anal.prior[i]);
  //  }

  //  printf("finished with priors\n");
  //

  // parms array declared
  //  int numParms = sizeof(pr)/sizeof(pr[0])/prCols;
  // double parms[numParms];
  double *parms = new double[numParms];

  // declare analysis
  anal.Y.assign(Y, Y + numDataRows);
  anal.n = numDataRows;
  if (suffStat) {
    anal.n_group.assign(N, N + numDataRows);
    anal.sd.assign(SD, SD + numDataRows);
  }
  anal.doses.assign(D, D + numDataRows);
  anal.disttype = dist;
  if (!detectAdvDir) {
    anal.isIncreasing = isIncreasing;
  }

  anal.alpha = alpha;
  anal.BMD_type = BMD_type;  // 1=absdev, 2 = stddev, 3 = reldev, 4 = pt, 5 = extra, 6 =
                             // hybrid_extra, 7 = hybrid_added   from src/include/cmodeldefs.h
  anal.BMR = BMRF;
  anal.samples = 0;  // num MCMC samples
  anal.tail_prob = 0.01;
  anal.suff_stat = suffStat;
  anal.parms = numParms;
  anal.prior_cols = prCols;
  anal.transform_dose = 0;
  anal.prior.assign(prior, prior + anal.prior_cols * anal.parms);
  anal.restricted = restricted;
  anal.detectAdvDir = detectAdvDir;

  printf("prior b4 adj:\n");
  for (int i = 0; i < prCols * numParms; i++) {
    printf("%.9f\n", anal.prior[i]);
  }

  struct python_continuous_model_result res;
  res.model = anal.model;
  res.nparms = anal.parms;
  res.dist_numE = 100;

  struct BMDS_results BMDSres;
  // set all parms as unbounded initially
  for (int i = 0; i < anal.parms; i++) {
    BMDSres.bounded.push_back(false);
    ;
    BMDSres.stdErr.push_back(BMDS_MISSING);
    BMDSres.lowerConf.push_back(BMDS_MISSING);
    BMDSres.upperConf.push_back(BMDS_MISSING);
  }
  BMDSres.BMD = -9999.0;
  BMDSres.BMDU = -9999.0;
  BMDSres.BMDL = -9999.0;
  BMDSres.AIC = -9999.0;
  res.bmdsRes = BMDSres;

  struct continuous_GOF gof;
  int nGOF;
  if (anal.suff_stat) {
    nGOF = anal.n;
  } else {
    // determine number of unique dose groups
    nGOF = 1;
    for (int i = 1; i < anal.n; i++) {
      int j = 0;
      for (j = 0; j < i; j++) {
        if (anal.doses[i] == anal.doses[j]) break;
      }
      if (i == j) nGOF++;
    }
  }
  res.gof = gof;

  struct continuous_AOD aod;
  std::vector<double> LL(5);
  std::vector<int> nParms(5);
  std::vector<double> AIC(5);
  double addConst;  // = 22.2;
  aod.LL = LL;
  aod.nParms = nParms;
  aod.AIC = AIC;
  aod.addConst = addConst;

  // double llRatio[4];
  // double DF[4];
  // double pVal[4];
  std::vector<double> llRatio(4);
  std::vector<double> DF(4);
  std::vector<double> pVal(4);
  struct testsOfInterest TOI;

  TOI.llRatio = llRatio;
  TOI.DF = DF;
  TOI.pVal = pVal;
  aod.TOI = TOI;
  res.aod = aod;

  printf("\n\n-----------INPUT---------------\n");
  printBmdsStruct(&anal);
  printf("\n\n");

  printf("calling pythonBMDSCont\n");
  pythonBMDSCont(&anal, &res);

  printf("\n\n");
  printf("prior after adj by model code:\n");
  for (int i = 0; i < prCols * numParms; i++) {
    printf("%.20f\n", anal.prior[i]);
  }

  printf("\n\n----------OUTPUT-----------\n");
  if (BMDSres.validResult || showResultsOverride) {
    printBmdsStruct(&res);
    //    printf("\nBMD Dist:\n");
    //    for (int i = 0; i < res.dist_numE; i++) {
    //      printf("i:%d, perc:%f, dist:%f\n", i, res.bmd_dist[i + res.dist_numE], res.bmd_dist[i]);
    //    }
    //
  } else {
    printf("Model was not run\n");
  }

  // debugContAnal(&anal);
}

void runPythonMultitumorAnalysis() {
  enum dich_model model = d_multistage;  // d_hill =1, d_gamma=2,d_logistic=3, d_loglogistic=4,
                                         // d_logprobit=5, d_multistage=6,d_probit=7,
                                         // d_qlinear=8,d_weibull=9
  int modelType = 1;                     // 1 = frequentist, 2 = bayesian
  bool restricted = true;                // only used for frequentist models
  int BMD_type = 1;                      // 1 = extra ; added otherwise
  double BMR = 0.1;
  double alpha = 0.05;
  std::vector<std::vector<double>> doses;
  std::vector<std::vector<double>> Y;
  std::vector<std::vector<double>> n_group;

  // data 1st test
  //  std::vector<double> doses1 = {0,50,100,150,200};
  //  std::vector<double> Y1 = {0,5,30,65,90};
  //  std::vector<double> n_group1 = {100,100,100,100,100};
  //  std::vector<double> doses2 = {0,50,100,150,200};
  //  std::vector<double> Y2 = {5,10,33,67,93};
  //  std::vector<double> n_group2 = {100,100,100,100,100};
  //  std::vector<double> doses3 = {0,50,100,150,200};
  //  std::vector<double> Y3 = {1,68,78,88,98};
  //  std::vector<double> n_group3 = {100,100,100,100,100};
  //  std::vector<int> n = {5,5,5};
  //  std::vector<int> degree = {0,0,0};
  //  doses.push_back(doses1);
  //  doses.push_back(doses2);
  //  doses.push_back(doses3);
  //  Y.push_back(Y1);
  //  Y.push_back(Y2);
  //  Y.push_back(Y3);
  //  n_group.push_back(n_group1);
  //  n_group.push_back(n_group2);
  //  n_group.push_back(n_group3);

  //  // online test
  //  std::vector<double> doses1 = {0, 50, 100, 200, 400};
  //  std::vector<double> Y1 = {0, 1, 2, 10, 19};
  //  std::vector<double> n_group1 = {20, 20, 20, 20, 20};
  //  std::vector<double> doses2 = {0, 50, 100, 200, 400};
  //  std::vector<double> Y2 = {0, 1, 2, 4, 11};
  //  std::vector<double> n_group2 = {20, 20, 20, 20, 20};
  //  std::vector<double> doses3 = {0, 50, 100, 200, 400};
  //  std::vector<double> Y3 = {0, 2, 2, 6, 9};
  //  std::vector<double> n_group3 = {20, 20, 20, 20, 20};
  //  std::vector<int> n = {5, 5, 5};
  //  std::vector<int> degree = {0, 0, 0};
  //  // std::vector<int> degree = {2,2,2};
  //  doses.push_back(doses1);
  //  doses.push_back(doses2);
  //  doses.push_back(doses3);
  //  Y.push_back(Y1);
  //  Y.push_back(Y2);
  //  Y.push_back(Y3);
  //  n_group.push_back(n_group1);
  //  n_group.push_back(n_group2);
  //  n_group.push_back(n_group3);

  // data test one recommended model
  //  std::vector<double> doses1 = {0,25,75,125,200};
  //  std::vector<double> Y1 = {0,0,1,7,11};
  //  std::vector<double> n_group1 = {20,20,20,20,20};
  //  std::vector<double> doses2 = {0,50,100,200,400};
  //  std::vector<double> Y2 = {1,68,78,88,98};
  //  std::vector<double> n_group2 = {100,100,100,100,100};
  //  std::vector<int> n = {5,5};
  //  std::vector<int> degree = {0,0};
  //  doses.push_back(doses1);
  //  doses.push_back(doses2);
  //  Y.push_back(Y1);
  //  Y.push_back(Y2);
  //  n_group.push_back(n_group1);
  //  n_group.push_back(n_group2);

  // data test one dataset
  //  std::vector<double> doses1 = {0,50,100,150,200};
  //  std::vector<double> Y1 = {0,5,30,65,90};
  //  std::vector<double> n_group1 = {100,100,100,100,100};
  //  std::vector<int> n = {5};
  //  std::vector<int> degree = {0};
  //  doses.push_back(doses1);
  //  Y.push_back(Y1);
  //  n_group.push_back(n_group1);

  // data test one recommended model
  //  std::vector<double> doses1 = {0,50,100,150,200};
  //  std::vector<double> Y1 = {0,5,30,65,90};
  //  std::vector<double> n_group1 = {100,100,100,100,100};
  //  std::vector<double> doses2 = {0,50,100,200,400};
  //  std::vector<double> Y2 = {1,68,78,88,98};
  //  std::vector<double> n_group2 = {100,100,100,100,100};
  //  std::vector<int> n = {5,5};
  //  std::vector<int> degree = {0,0};
  //  doses.push_back(doses1);
  //  doses.push_back(doses2);
  //  Y.push_back(Y1);
  //  Y.push_back(Y2);
  //  n_group.push_back(n_group1);
  //  n_group.push_back(n_group2);

  //  //data test no recommended model
  //  std::vector<double> doses1 = {0,50,100,200,400};
  //  std::vector<double> Y1 = {1,68,78,88,98};
  //  std::vector<double> n_group1 = {100,100,100,100,100};
  //  std::vector<int> n = {5};
  //  std::vector<int> degree = {0};
  //  doses.push_back(doses1);
  //  Y.push_back(Y1);
  //  n_group.push_back(n_group1);

  // MT test dataset #1 2 well fitting datasets
  std::vector<double> doses1 = {0, 50, 100, 200, 400};
  std::vector<double> Y1 = {0, 1, 4, 8, 16};
  std::vector<double> n_group1 = {50, 50, 50, 50, 50};
  std::vector<double> doses2 = {0, 50, 100, 200, 400};
  std::vector<double> Y2 = {0, 4, 12, 36, 48};
  std::vector<double> n_group2 = {50, 50, 50, 50, 50};
  std::vector<int> n = {5, 5};
  std::vector<int> degree = {0, 0};
  doses.push_back(doses1);
  doses.push_back(doses2);
  Y.push_back(Y1);
  Y.push_back(Y2);
  n_group.push_back(n_group1);
  n_group.push_back(n_group2);

  //  // MT test dataset #2  3 well-fitting datasets
  //  std::vector<double> doses1 = {0, 50, 100, 150, 200};
  //  std::vector<double> Y1 = {0, 0, 4, 8, 16};
  //  std::vector<double> n_group1 = {50, 50, 50, 50, 50};
  //  std::vector<double> doses2 = {0, 50, 100, 150, 200};
  //  std::vector<double> Y2 = {0, 4, 12, 36, 48};
  //  std::vector<double> n_group2 = {50, 50, 50, 50, 50};
  //  std::vector<double> doses3 = {0, 50, 100, 150, 200};
  //  std::vector<double> Y3 = {0, 2, 3, 5, 8};
  //  std::vector<double> n_group3 = {50, 50, 50, 50, 50};
  //  std::vector<int> n = {5,5,5};
  //  std::vector<int> degree = {0, 0, 0};
  //  doses.push_back(doses1);
  //  doses.push_back(doses2);
  //  doses.push_back(doses3);
  //  Y.push_back(Y1);
  //  Y.push_back(Y2);
  //  Y.push_back(Y3);
  //  n_group.push_back(n_group1);
  //  n_group.push_back(n_group2);
  //  n_group.push_back(n_group3);

  //  // MT test dataset #3 2 well-fitting and 1 poorly fitting datasets
  //  std::vector<double> doses1 = {0, 50, 100, 150, 200};
  //  std::vector<double> Y1 = {0, 0, 4, 8, 16};
  //  std::vector<double> n_group1 = {50, 50, 50, 50, 50};
  //  std::vector<double> doses2 = {0, 50, 100, 150, 200};
  //  std::vector<double> Y2 = {0, 3, 12, 20, 36};
  //  std::vector<double> n_group2 = {50, 50, 50, 50, 50};
  //  std::vector<double> doses3 = {0, 50, 100, 150, 200};
  //  std::vector<double> Y3 = {0, 2, 3, 5, 8};
  //  std::vector<double> n_group3 = {50, 50, 50, 50, 50};
  //  std::vector<int> n = {5,5,5};
  //  std::vector<int> degree = {0, 0, 0};
  //  doses.push_back(doses1);
  //  doses.push_back(doses2);
  //  doses.push_back(doses3);
  //  Y.push_back(Y1);
  //  Y.push_back(Y2);
  //  Y.push_back(Y3);
  //  n_group.push_back(n_group1);
  //  n_group.push_back(n_group2);
  //  n_group.push_back(n_group3);

  //  // MT test dataset #4  3 well-fitting datasets, 1 with high dose dropped
  //  std::vector<double> doses1 = {0, 50, 100, 150, 200};
  //  std::vector<double> Y1 = {0, 0, 4, 8, 16};
  //  std::vector<double> n_group1 = {50, 50, 50, 50, 50};
  //  std::vector<double> doses2 = {0, 50, 100, 150};
  //  std::vector<double> Y2 = {0, 3, 12, 20};
  //  std::vector<double> n_group2 = {50, 50, 50, 50};
  //  std::vector<double> doses3 = {0, 50, 100, 150, 200};
  //  std::vector<double> Y3 = {0, 2, 3, 5, 8};
  //  std::vector<double> n_group3 = {50, 50, 50, 50, 50};
  //  std::vector<int> n = {5,4,5};
  //  std::vector<int> degree = {0, 0, 0};
  //  doses.push_back(doses1);
  //  doses.push_back(doses2);
  //  doses.push_back(doses3);
  //  Y.push_back(Y1);
  //  Y.push_back(Y2);
  //  Y.push_back(Y3);
  //  n_group.push_back(n_group1);
  //  n_group.push_back(n_group2);
  //  n_group.push_back(n_group3);

  //  //MT test dataset #5 2 well fitting datasets, 2 poorly fitting datasets
  //  std::vector<double> doses1 = {0, 100, 200, 400, 800};
  //  std::vector<double> Y1 = {0,3,17,41,49};
  //  std::vector<double> n_group1 = {50,50,50,50,50};
  //  std::vector<double> doses2 = {0, 100, 200, 400, 800};
  //  std::vector<double> Y2 = {0,3,6,12,24};
  //  std::vector<double> n_group2 = {50,50,50,50,50};
  //  std::vector<double> doses3 = {0, 100, 200, 400, 800};
  //  std::vector<double> Y3 = {0,3,6,12,8};
  //  std::vector<double> n_group3 = {50,50,50,50,50};
  //  std::vector<double> doses4 = {0, 100, 200, 400, 800};
  //  std::vector<double> Y4 = {0,1,2,4,49};
  //  std::vector<double> n_group4 = {50,50,50,50,50};
  //  std::vector<int> n = {5,5,5,5};
  //  std::vector<int> degree = {0, 0, 0, 0};
  //  doses.push_back(doses1);
  //  doses.push_back(doses2);
  //  doses.push_back(doses3);
  //  doses.push_back(doses4);
  //  Y.push_back(Y1);
  //  Y.push_back(Y2);
  //  Y.push_back(Y3);
  //  Y.push_back(Y4);
  //  n_group.push_back(n_group1);
  //  n_group.push_back(n_group2);
  //  n_group.push_back(n_group3);
  //  n_group.push_back(n_group4);

  //  //MT test dataset #6 4 well fitting datasets, 2 with high responses and 2 with low responses
  //  std::vector<double> doses1 = {0, 100, 200, 400, 800};
  //  std::vector<double> Y1 = {0,5,15,30,45};
  //  std::vector<double> n_group1 = {50,50,50,50,50};
  //  std::vector<double> doses2 = {0, 100, 200, 400, 800};
  //  std::vector<double> Y2 = {0,3,23,49,50};
  //  std::vector<double> n_group2 = {50,50,50,50,50};
  //  std::vector<double> doses3 = {0, 100, 200, 400, 800};
  //  std::vector<double> Y3 = {0,0,2,3,6};
  //  std::vector<double> n_group3 = {50,50,50,50,50};
  //  std::vector<double> doses4 = {0, 100, 200, 400, 800};
  //  std::vector<double> Y4 = {0,1,2,4,8};
  //  std::vector<double> n_group4 = {50,50,50,50,50};
  //  std::vector<int> n = {5,5,5,5};
  //  std::vector<int> degree = {0, 0, 0, 0};
  //  doses.push_back(doses1);
  //  doses.push_back(doses2);
  //  doses.push_back(doses3);
  //  doses.push_back(doses4);
  //  Y.push_back(Y1);
  //  Y.push_back(Y2);
  //  Y.push_back(Y3);
  //  Y.push_back(Y4);
  //  n_group.push_back(n_group1);
  //  n_group.push_back(n_group2);
  //  n_group.push_back(n_group3);
  //  n_group.push_back(n_group4);

  //  //MT test dataset #7 4 well fitting datasets with different Ns
  //  std::vector<double> doses1 = {0, 100, 200, 400, 800};
  //  std::vector<double> Y1 = {0,1,2,5,12};
  //  std::vector<double> n_group1 = {50,50,50,50,50};
  //  std::vector<double> doses2 = {0, 100, 200, 400, 800};
  //  std::vector<double> Y2 = {0,1,3,6,18};
  //  std::vector<double> n_group2 = {50,50,50,50,50};
  //  std::vector<double> doses3 = {0, 100, 200, 400, 800};
  //  std::vector<double> Y3 = {0,0,1,5,12};
  //  std::vector<double> n_group3 = {25,25,25,25,25};
  //  std::vector<double> doses4 = {0, 100, 200, 400, 800};
  //  std::vector<double> Y4 = {0,1,4,21,24};
  //  std::vector<double> n_group4 = {25,25,25,25,25};
  //  std::vector<int> n = {5,5,5,5};
  //  std::vector<int> degree = {0, 0, 0, 0};
  //  doses.push_back(doses1);
  //  doses.push_back(doses2);
  //  doses.push_back(doses3);
  //  doses.push_back(doses4);
  //  Y.push_back(Y1);
  //  Y.push_back(Y2);
  //  Y.push_back(Y3);
  //  Y.push_back(Y4);
  //  n_group.push_back(n_group1);
  //  n_group.push_back(n_group2);
  //  n_group.push_back(n_group3);
  //  n_group.push_back(n_group4);

  //  //MT test dataset #8 3 well fitting datasets with different Ns and 1 with high background
  //  response std::vector<double> doses1 = {0, 100, 200, 400, 800}; std::vector<double> Y1 =
  //  {0,1,2,9,23}; std::vector<double> n_group1 = {50,50,50,50,50}; std::vector<double> doses2 =
  //  {0, 100, 200, 400, 800}; std::vector<double> Y2 = {5,1,3,6,18}; std::vector<double> n_group2 =
  //  {50,50,50,50,50}; std::vector<double> doses3 = {0, 100, 200, 400, 800}; std::vector<double> Y3
  //  = {5,7,12,18,25}; std::vector<double> n_group3 = {25,25,25,25,25}; std::vector<int> n =
  //  {5,5,5}; std::vector<int> degree = {0, 0, 0}; doses.push_back(doses1);
  //  doses.push_back(doses2);
  //  doses.push_back(doses3);
  //  Y.push_back(Y1);
  //  Y.push_back(Y2);
  //  Y.push_back(Y3);
  //  n_group.push_back(n_group1);
  //  n_group.push_back(n_group2);
  //  n_group.push_back(n_group3);

  //  //MT test dataset #9 3 well fitting datasets, 1 poorly fitting dataset, 1 well fitting with
  //  high dose dropped std::vector<double> doses1 = {0, 100, 200, 400, 800}; std::vector<double> Y1
  //  = {0,1,2,9,22}; std::vector<double> n_group1 = {50,50,50,50,50}; std::vector<double> doses2 =
  //  {0, 100, 200, 400, 800}; std::vector<double> Y2 = {1,1,3,6,18}; std::vector<double> n_group2 =
  //  {50,50,50,50,50}; std::vector<double> doses3 = {0, 100, 200, 400, 800}; std::vector<double> Y3
  //  = {0,1,5,11,19}; std::vector<double> n_group3 = {50,50,50,50,50}; std::vector<double> doses4 =
  //  {0, 100, 200, 400, 800}; std::vector<double> Y4 = {0,1,11,2,22}; std::vector<double> n_group4
  //  = {50,50,50,50,50}; std::vector<double> doses5 = {0, 100, 200, 400}; std::vector<double> Y5 =
  //  {0,1,6,14}; std::vector<double> n_group5 = {50,50,50,50}; std::vector<int> n = {5,5,5,5,4};
  //  std::vector<int> degree = {0, 0, 0, 0, 0};
  //  doses.push_back(doses1);
  //  doses.push_back(doses2);
  //  doses.push_back(doses3);
  //  doses.push_back(doses4);
  //  doses.push_back(doses5);
  //  Y.push_back(Y1);
  //  Y.push_back(Y2);
  //  Y.push_back(Y3);
  //  Y.push_back(Y4);
  //  Y.push_back(Y5);
  //  n_group.push_back(n_group1);
  //  n_group.push_back(n_group2);
  //  n_group.push_back(n_group3);
  //  n_group.push_back(n_group4);
  //  n_group.push_back(n_group5);

  //  //MT test dataset #10 8 well fitting datasets, 1 well fitting with high dose dropped, 1 well
  //  fitting with different doses std::vector<double> doses1 = {0, 100, 200, 400, 800};
  //  std::vector<double> Y1 = {0,1,2,9,22};
  //  std::vector<double> n_group1 = {50,50,50,50,50};
  //  std::vector<double> doses2 = {0, 100, 200, 400, 800};
  //  std::vector<double> Y2 = {1,1,3,6,18};
  //  std::vector<double> n_group2 = {50,50,50,50,50};
  //  std::vector<double> doses3 = {0, 100, 200, 400, 800};
  //  std::vector<double> Y3 = {0,1,5,11,19};
  //  std::vector<double> n_group3 = {50,50,50,50,50};
  //  std::vector<double> doses4 = {0, 50, 100, 200, 500};
  //  std::vector<double> Y4 = {0,1,3,9,22};
  //  std::vector<double> n_group4 = {50,50,50,50,50};
  //  std::vector<double> doses5 = {0, 100, 200, 400};
  //  std::vector<double> Y5 = {0,1,6,14};
  //  std::vector<double> n_group5 = {50,50,50,50};
  //  std::vector<double> doses6 = {0, 100, 200, 400, 800};
  //  std::vector<double> Y6 = {1,4,8,19,40};
  //  std::vector<double> n_group6 = {50,50,50,50,50};
  //  std::vector<double> doses7 = {0, 100, 200, 400, 800};
  //  std::vector<double> Y7 = {1,1,3,4,21};
  //  std::vector<double> n_group7 = {50,50,50,50,50};
  //  std::vector<double> doses8 = {0, 100, 200, 400, 800};
  //  std::vector<double> Y8 = {0,1,4,21,45};
  //  std::vector<double> n_group8 = {50,50,50,50,50};
  //  std::vector<double> doses9 = {0, 100, 200, 400, 800};
  //  std::vector<double> Y9 = {0,1,8,32,45};
  //  std::vector<double> n_group9 = {50,50,50,50,50};
  //  std::vector<double> doses10 = {0, 100, 200, 400, 800};
  //  std::vector<double> Y10 = {0,0,5,29,49};
  //  std::vector<double> n_group10 = {50,50,50,50,50};
  //  std::vector<int> n = {5,5,5,5,4,5,5,5,5,5,5};
  //  std::vector<int> degree = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  //  doses.push_back(doses1);
  //  doses.push_back(doses2);
  //  doses.push_back(doses3);
  //  doses.push_back(doses4);
  //  doses.push_back(doses5);
  //  doses.push_back(doses6);
  //  doses.push_back(doses7);
  //  doses.push_back(doses8);
  //  doses.push_back(doses9);
  //  doses.push_back(doses10);
  //  Y.push_back(Y1);
  //  Y.push_back(Y2);
  //  Y.push_back(Y3);
  //  Y.push_back(Y4);
  //  Y.push_back(Y5);
  //  Y.push_back(Y6);
  //  Y.push_back(Y7);
  //  Y.push_back(Y8);
  //  Y.push_back(Y9);
  //  Y.push_back(Y10);
  //  n_group.push_back(n_group1);
  //  n_group.push_back(n_group2);
  //  n_group.push_back(n_group3);
  //  n_group.push_back(n_group4);
  //  n_group.push_back(n_group5);
  //  n_group.push_back(n_group6);
  //  n_group.push_back(n_group7);
  //  n_group.push_back(n_group8);
  //  n_group.push_back(n_group9);
  //  n_group.push_back(n_group10);

  /////////////////////////////////////////////////
  ////END USER INPUT
  ////////////////////////////////////////////////////

  // priors defined columnwise
  int prCols = 5;

  int numDatasets = Y.size();
  // int numDatasets = 1;

  struct python_multitumor_analysis anal;
  anal.ndatasets = numDatasets;
  anal.n = n;
  anal.degree = degree;
  anal.BMR = BMR;
  anal.BMD_type = BMD_type;
  anal.alpha = alpha;

  struct python_multitumor_result res;
  res.ndatasets = numDatasets;

  // create individual models analyses
  std::vector<std::vector<python_dichotomous_analysis>> models;
  int count = 0;

  for (int dataset = 0; dataset < numDatasets; dataset++) {
    int numDataRows = Y[dataset].size();

    struct python_dichotomous_analysis modAnal;
    modAnal.model = d_multistage;
    printf("model = %d\n", modAnal.model);
    modAnal.BMD_type = BMD_type;
    modAnal.BMR = BMR;
    modAnal.alpha = alpha;
    modAnal.prior_cols = prCols;

    struct python_dichotomous_model_result modRes;
    modRes.model = modAnal.model;
    modRes.dist_numE = 200;

    std::vector<python_dichotomous_analysis> modGroup;
    std::vector<python_dichotomous_model_result> modResGroup;

    // Individual model construction
    // declare analysis
    modAnal.parms = 1 + degree[dataset];
    modAnal.Y = Y[dataset];
    modAnal.n_group = n_group[dataset];
    modAnal.doses = doses[dataset];
    modAnal.n = numDataRows;

    // needs to be changed based on model degree
    if (degree[dataset] == 0) {
      // handle autoselect degree
      count = 0;
      // run models from degree = 1 to k-1, where k = number of dose groups
      for (int deg = 1; deg < anal.n[dataset]; deg++) {
        modAnal.degree = deg;
        modAnal.prior = getMultitumorPrior(modAnal.degree, modAnal.prior_cols);
        modAnal.parms = modAnal.degree + 1;
        modRes.nparms = modAnal.parms;
        modGroup.push_back(modAnal);
        modResGroup.push_back(modRes);
        count++;
      }
    } else {
      modAnal.degree = degree[dataset];
      modAnal.prior = getMultitumorPrior(modAnal.degree, modAnal.prior_cols);
      modRes.nparms = modAnal.parms;
      modGroup.push_back(modAnal);
      modResGroup.push_back(modRes);
      count = 1;
    }

    anal.nmodels.push_back(count);
    anal.models.push_back(modGroup);
    res.nmodels.push_back(count);
    res.models.push_back(modResGroup);
  }

  // run MSCombo
  pythonBMDSMultitumor(&anal, &res);

  std::cout << std::endl;
  if (res.validResult) {
    std::cout << "INPUT" << std::endl;
    printBmdsStruct(&anal);
    std::cout << "OUTPUT" << std::endl;
    printBmdsStruct(&res);
  } else {
    std::cout << "Multitumor analysis failed" << std::endl;
  }
}

void runTestMultitumorModel() {
  enum dich_model model = d_multistage;  // d_hill =1, d_gamma=2,d_logistic=3, d_loglogistic=4,
                                         // d_logprobit=5, d_multistage=6,d_probit=7,
                                         // d_qlinear=8,d_weibull=9
  int modelType = 1;                     // 1 = frequentist, 2 = bayesian
  bool restricted = true;                // only used for frequentist models
  int BMD_type = 1;                      // 1 = extra ; added otherwise
  double BMR = 0.1;
  double alpha = 0.05;

  // data
  std::vector<double> doses1 = {0, 50, 100, 150, 200};
  std::vector<double> Y1 = {0, 5, 30, 65, 90};
  std::vector<double> n_group1 = {100, 100, 100, 100, 100};
  std::vector<double> doses2 = {0, 50, 100, 150, 200};
  std::vector<double> Y2 = {5, 10, 33, 67, 93};
  std::vector<double> n_group2 = {100, 100, 100, 100, 100};
  std::vector<double> doses3 = {0, 50, 100, 150, 200};
  std::vector<double> Y3 = {1, 68, 78, 88, 98};
  std::vector<double> n_group3 = {100, 100, 100, 100, 100};

  std::vector<std::vector<double>> doses;
  std::vector<std::vector<double>> Y;
  std::vector<std::vector<double>> n_group;
  doses.push_back(doses1);
  doses.push_back(doses2);
  doses.push_back(doses3);
  Y.push_back(Y1);
  Y.push_back(Y2);
  Y.push_back(Y3);
  n_group.push_back(n_group1);
  n_group.push_back(n_group2);
  n_group.push_back(n_group3);

  std::vector<int> n = {5, 5, 5};
  std::vector<int> degree = {2, 2, 2};
  /////////////////////////////////////////////////
  ////END USER INPUT
  ////////////////////////////////////////////////////

  // priors defined columnwise
  int prCols = 5;

  int numDatasets = Y.size();
  // int numDatasets = 1;

  struct python_multitumor_analysis anal;
  anal.ndatasets = numDatasets;
  anal.n = n;
  anal.degree = degree;
  anal.BMR = BMR;
  anal.BMD_type = BMD_type;
  anal.alpha = alpha;

  struct python_multitumor_result res;
  res.ndatasets = numDatasets;

  // create individual models analyses
  std::vector<std::vector<python_dichotomous_analysis>> models;
  int count = 0;

  for (int dataset = 0; dataset < numDatasets; dataset++) {
    std::cout << "adding dataset:" << dataset << std::endl;

    int numDataRows = Y[dataset].size();

    struct python_dichotomous_analysis modAnal;
    modAnal.model = d_multistage;
    printf("model = %d\n", modAnal.model);
    modAnal.BMD_type = BMD_type;
    modAnal.BMR = BMR;
    modAnal.alpha = alpha;
    modAnal.prior_cols = prCols;

    struct python_dichotomous_model_result modRes;
    modRes.model = modAnal.model;
    modRes.dist_numE = 200;

    std::vector<python_dichotomous_analysis> modGroup;
    std::vector<python_dichotomous_model_result> modResGroup;

    // Individual model construction
    // declare analysis
    modAnal.parms = 1 + degree[dataset];
    modAnal.Y = Y[dataset];
    modAnal.n_group = n_group[dataset];
    modAnal.doses = doses[dataset];
    modAnal.n = numDataRows;

    // needs to be changed based on model degree
    if (degree[dataset] == 0) {
      // handle autoselect degree
      count = 0;
      for (int deg = 2; deg < anal.n[dataset]; deg++) {
        modAnal.degree = deg;
        modAnal.prior = getMultitumorPrior(modAnal.degree, modAnal.prior_cols);
        modAnal.parms = modAnal.degree + 1;
        modRes.nparms = modAnal.parms;
        modGroup.push_back(modAnal);
        modResGroup.push_back(modRes);
        count++;
      }
    } else {
      modAnal.degree = degree[dataset];
      modAnal.prior = getMultitumorPrior(modAnal.degree, modAnal.prior_cols);
      modRes.nparms = modAnal.parms;
      modGroup.push_back(modAnal);
      modResGroup.push_back(modRes);
      count = 1;
    }

    anal.nmodels.push_back(count);
    anal.models.push_back(modGroup);
    res.nmodels.push_back(count);
    res.models.push_back(modResGroup);
  }

  //  std::cout<<"adding selected model"<<std::endl;
  //  res.selectedModelIndex.push_back(0);
  //  res.selectedModelIndex.push_back(0);
  //  res.selectedModelIndex.push_back(0);
  //  std::cout<<"Selected model Indexes:  ";

  //  res.validResult.push_back(true);
  //  res.validResult.push_back(true);
  //  res.validResult.push_back(true);

  //  std::cout<<"adding res[0]"<<std::endl;
  //  res.models[0][0].max = -182.8268;
  //  res.models[0][0].parms.push_back(0.0);
  //  res.models[0][0].parms.push_back(0.0);
  //  res.models[0][0].parms.push_back(4.56977e-05);
  //  res.models[0][0].bmd = 48.0167;
  //  res.models[0][0].bmdsRes.BMD = 48.0167;
  //  res.models[0][0].bmdsRes.BMDL = 44.1401;
  //  res.models[0][0].bmdsRes.BMDU = 51.2664;
  //  res.models[0][0].bmdsRes.AIC = 367.736;
  //  res.models[0][0].bmdsRes.chisq = 8.17;
  //
  //  std::cout<<"adding res[1]"<<std::endl;
  //  res.models[1][0].max = -209.332;
  //  res.models[1][0].parms.push_back(0.034571);
  //  res.models[1][0].parms.push_back(0.0);
  //  res.models[1][0].parms.push_back(4.90516e-05);
  //  res.models[1][0].bmd = 47.2147;
  //  res.models[1][0].bmdsRes.BMD = 47.2147;
  //  res.models[1][0].bmdsRes.BMDL = 42.6674;
  //  res.models[1][0].bmdsRes.BMDU = 50.8673;
  //  res.models[1][0].bmdsRes.AIC = 422.664;
  //  res.models[1][0].bmdsRes.chisq = 8.71;
  //
  //  std::cout<<"adding res[2]"<<std::endl;
  //  res.models[2][0].max = -171.622;
  //  res.models[2][0].parms.push_back(0.011097);
  //  res.models[2][0].parms.push_back(0.0169765);
  //  res.models[2][0].parms.push_back(0.0);
  //  res.models[2][0].bmd = 6.27976;
  //  res.models[2][0].bmdsRes.BMD = 6.27976;
  //  res.models[2][0].bmdsRes.BMDL = 5.59913;
  //  res.models[2][0].bmdsRes.BMDU = 7.40053;
  //  res.models[2][0].bmdsRes.AIC = 367.736;
  //  res.models[2][0].bmdsRes.chisq = 8.17;

  std::cout << "calling runMultitumor" << std::endl;
  // run MSCombo
  //  pythonBMDSMultitumor(&anal, &res);

  std::cout << "here" << std::endl;
  //  runMultitumorModel(&anal, &res);
  pythonBMDSMultitumor(&anal, &res);
  //  testCall();

  //  //individual model results
  //  for (int dataset=0; dataset<numDatasets; dataset++){
  //    std::cout<<"dataset:"<<dataset<<std::endl;
  //    for (int mod=0; mod<anal.nmodels[dataset]; mod++){
  //        std::cout<<" model:"<<mod<<std::endl;
  //        printDichoModResult(&anal.models[dataset][mod], &res.models[dataset][mod],true);
  //    }
  //  }

  std::cout << "BMD:  " << res.BMD << std::endl;
  std::cout << "BMDL: " << res.BMDL << std::endl;
  std::cout << "BMDU: " << res.BMDU << std::endl;
  std::cout << "slope factor: " << res.slopeFactor << std::endl;
  std::cout << "combined LL: " << res.combined_LL << std::endl;
  std::cout << "combined LL constant: " << res.combined_LL_const << std::endl;
}

void printDichoModResult(
    struct python_dichotomous_analysis *pyAnal, struct python_dichotomous_model_result *pyRes,
    bool showResultsOverride
) {
  printf("tlink bmdsRes.validResult = %s\n", pyRes->bmdsRes.validResult ? "valid" : "invalid");
  if (pyRes->bmdsRes.validResult || showResultsOverride) {
    std::cout << "Valid Result" << std::endl;
    printf("\nBenchmark Dose\n");
    printf("max: %f\n", pyRes->max);
    printf("BMD: %f\n", pyRes->bmdsRes.BMD);
    printf("BMDL: %f\n", pyRes->bmdsRes.BMDL);
    printf("BMDU: %f\n", pyRes->bmdsRes.BMDU);
    printf("AIC: %f\n", pyRes->bmdsRes.AIC);
    printf("LPP: %f\n", pyRes->bmdsRes.BIC_equiv);
    printf("P-value: %f\n", pyRes->gof.p_value);
    printf("DOF: %f\n", pyRes->gof.df);
    printf("Chi^2: %f\n", pyRes->bmdsRes.chisq);

    printf("\nModel Parameters\n");
    printf("# of parms: %d\n", pyAnal->parms);
    printf("parm, estimate, bounded, std.err., lower conf, upper conf\n");
    for (int i = 0; i < pyAnal->parms; i++) {
      printf(
          "%d, %.10f, %s, %f, %f, %f\n", i, pyRes->parms[i],
          pyRes->bmdsRes.bounded[i] ? "true" : "false", pyRes->bmdsRes.stdErr[i],
          pyRes->bmdsRes.lowerConf[i], pyRes->bmdsRes.upperConf[i]
      );
    }

    printf("\ncov matrix\n");
    for (int i = 0; i < pyAnal->parms * pyAnal->parms; i++) {
      printf("%d, %f\n", i, pyRes->cov[i]);
    }

    printf("\nGoodness of Fit\n");
    printf("Dose, EstProb, Expected, Observed, Size, ScaledRes\n");
    for (int i = 0; i < pyRes->gof.n; i++) {
      printf(
          "%f, %f, %f, %f, %f, %f\n", pyAnal->doses[i], pyRes->gof.expected[i] / pyAnal->n_group[i],
          pyRes->gof.expected[i], pyAnal->Y[i], pyAnal->n_group[i], pyRes->gof.residual[i]
      );
    }
    printf("\nError bars\n");
    for (int i = 0; i < pyRes->gof.n; i++) {
      printf("%f, %f\n", pyRes->gof.ebLower[i], pyRes->gof.ebUpper[i]);
    }

    printf("\nAnalysis of Deviance\n");
    printf("  Model,   LL,    #parms,   deviance,   test DF,  pval\n");
    printf("Full Model,  %f,  %d,  -,  -,  NA\n", pyRes->aod.fullLL, pyRes->aod.nFull);
    printf(
        "Fitted Model,  %f,  %d,  %f,  %d,  %f\n", pyRes->aod.fittedLL, pyRes->aod.nFit,
        pyRes->aod.devFit, pyRes->aod.dfFit, pyRes->aod.pvFit
    );
    printf(
        "Reduced Model,  %f,  %d,  %f,  %d,  %f\n", pyRes->aod.redLL, pyRes->aod.nRed,
        pyRes->aod.devRed, pyRes->aod.dfRed, pyRes->aod.pvRed
    );

    printf("\nBMD Dist:\n");
    for (int i = 0; i < pyRes->dist_numE; i++) {
      printf(
          "i:%d, perc:%f, dist:%f\n", i, pyRes->bmd_dist[i + pyRes->dist_numE], pyRes->bmd_dist[i]
      );
    }
  } else {
    printf("\nModel was not run\n");
  }
}
std::vector<double> getMultitumorPrior(int degree, int prior_cols) {
  int numParms = 2;  // 2 parameters for multistage cancer G & B
  // initial values for multistage 2
  std::vector<double> prG(prRFreqMultistageCancerG, prRFreqMultistageCancerG + prior_cols);
  std::vector<double> prB(prRFreqMultistageCancerB, prRFreqMultistageCancerB + prior_cols);

  std::vector<double> pr;
  for (int i = 0; i < prior_cols; i++) {
    pr.push_back(prG[i]);
    for (int j = 0; j < degree; j++) {
      pr.push_back(prB[i]);
    }
  }
  return pr;
}

std::vector<double> getNlogisticPrior(int ngrp, int prior_cols, bool restricted) {
  std::vector<double> prG(prNLogisticG, prNLogisticG + prior_cols);
  std::vector<double> prB(prNLogisticB, prNLogisticB + prior_cols);
  std::vector<double> prT1(prNLogisticT1, prNLogisticT1 + prior_cols);
  std::vector<double> prT2(prNLogisticT2, prNLogisticT2 + prior_cols);
  std::vector<double> prURho(prUNLogisticRho, prUNLogisticRho + prior_cols);
  std::vector<double> prRRho(prRNLogisticRho, prRNLogisticRho + prior_cols);
  std::vector<double> prPhi(prNLogisticPhi, prNLogisticPhi + prior_cols);

  std::vector<double> pr;

  for (int i = 0; i < prior_cols; i++) {
    pr.push_back(prG[i]);
    pr.push_back(prB[i]);
    pr.push_back(prT1[i]);
    pr.push_back(prT2[i]);
    if (restricted) {
      pr.push_back(prRRho[i]);
    } else {
      pr.push_back(prURho[i]);
    }
    for (int j = 0; j < ngrp; j++) {
      pr.push_back(prPhi[i]);
    }
  }

  return pr;
}

std::vector<double> getNctrPrior(int ngrp, int prior_cols, bool restricted) {
  std::vector<double> prG(prNctrG, prNctrG + prior_cols);
  std::vector<double> prB(prNctrB, prNctrB + prior_cols);

  // T1 and T2 get initialized based on smin/smax from data
  std::vector<double> prT1(prNctrT1, prNctrT1 + prior_cols);
  std::vector<double> prT2(prNctrT2, prNctrT2 + prior_cols);
  std::vector<double> prURho(prUNctrRho, prUNctrRho + prior_cols);
  std::vector<double> prRRho(prRNctrRho, prRNctrRho + prior_cols);
  std::vector<double> prPhi(prNctrPhi, prNctrPhi + prior_cols);

  std::vector<double> pr;

  for (int i = 0; i < prior_cols; i++) {
    pr.push_back(prG[i]);
    pr.push_back(prB[i]);
    pr.push_back(prT1[i]);
    pr.push_back(prT2[i]);
    //    pr.push_back(-9999);
    //    pr.push_back(-9999);
    if (restricted) {
      pr.push_back(prRRho[i]);
    } else {
      pr.push_back(prURho[i]);
    }
    for (int j = 0; j < ngrp; j++) {
      pr.push_back(prPhi[i]);
    }
  }

  return pr;
}

void runPythonNestedAnalysis() {
  bool showResultsOverride = true;
  struct python_nested_analysis pyAnal;
  // pyAnal.model = nlogistic;
  pyAnal.model = nctr;
  pyAnal.countAllParmsOnBoundary = true;
  bool isRestricted = true;

  // nested.dax
  //  std::vector<double> D = {0,  0,  0,  0,  0,   0,   0,   0,   0,   0,   25,  25,  25,
  //                           25, 25, 25, 25, 25,  25,  25,  50,  50,  50,  50,  50,  50,
  //                           50, 50, 50, 50, 100, 100, 100, 100, 100, 100, 100, 100, 100};
  //  std::vector<double> LS = {16, 9,  15, 14, 13, 9,  10, 14, 10, 11, 14, 9,  14,
  //                            9,  13, 12, 10, 10, 11, 14, 11, 11, 14, 11, 10, 11,
  //                            10, 15, 7,  14, 11, 14, 12, 13, 12, 14, 11, 8,  10};
  //  std::vector<double> I = {1, 1, 2, 3, 3, 0, 2, 2, 1, 2, 4, 5, 6, 2, 6, 3, 1, 2, 4, 3,
  //                           4, 5, 5, 4, 5, 4, 5, 6, 2, 4, 6, 6, 8, 7, 8, 6, 6, 5, 4};
  //  std::vector<double> LSC = {16, 9,  15, 14, 13, 9,  10, 14, 10, 11, 14, 9,  14,
  //                             9,  13, 12, 10, 10, 11, 14, 11, 11, 14, 11, 10, 11,
  //                             10, 15, 7,  14, 11, 14, 12, 13, 12, 14, 11, 8,  10};

  //  //  testing dataset id:1
  //  std::vector<double> D = {0,    0,    0,    0,    0,    0,    0,    0,   0,   0,   0,   0,   0,
  //                           0,    0,    0,    0,    0,    0,    0,    0,   30,  30,  30,  30, 30,
  //                           30,   30,   30,   30,   30,   30,   30,   30,  30,  30,  30,  30, 30,
  //                           30,   30,   300,  300,  300,  300,  300,  300, 300, 300, 300, 300,
  //                           300, 300,  300,  300,  300,  300,  300,  300,  300, 300, 300, 300,
  //                           300, 1200, 1200, 1200, 1200, 1200, 1200, 1200, 1200, 1200};
  //  std::vector<double> LS = {1,  3,  3,  5,  5,  5,  6,  7, 7, 8, 8, 8, 8, 8,  9,  9, 10, 10, 10,
  //                            10, 11, 3,  6,  6,  6,  7,  7, 8, 8, 8, 8, 9, 9,  9,  9, 9,  10, 10,
  //                            11, 11, 12, 2,  4,  5,  5,  7, 7, 7, 8, 8, 8, 9,  9,  9, 9,  10, 10,
  //                            10, 11, 11, 11, 11, 12, 16, 2, 3, 5, 7, 9, 9, 10, 10, 10};
  //  std::vector<double> I = {1, 0,  2, 0, 3, 1, 0,  0, 0, 0, 0, 1, 1, 2, 0, 0, 1, 0, 0,
  //                           0, 11, 2, 6, 4, 0, 4,  1, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 2,
  //                           0, 0,  0, 2, 4, 3, 5,  0, 1, 2, 8, 0, 0, 0, 1, 1, 1, 1, 1,
  //                           0, 0,  1, 0, 1, 0, 12, 2, 1, 0, 1, 1, 9, 2, 4, 1};
  //  std::vector<double> LSC = {1.0,  3.0,  3.0,  5.0,  5.0,  5.0,  6.0,  7.0,  7.0,  8.0,  8.0,
  //                             8.0,  8.0,  8.0,  9.0,  9.0,  10.0, 10.0, 10.0, 10.0, 11.0, 3.0,
  //                             6.0,  6.0,  6.0,  7.0,  7.0,  8.0,  8.0,  8.0,  8.0,  9.0,  9.0,
  //                             9.0,  9.0,  9.0,  10.0, 10.0, 11.0, 11.0, 12.0, 2.0,  4.0,  5.0,
  //                             5.0,  7.0,  7.0,  7.0,  8.0,  8.0,  8.0,  9.0,  9.0,  9.0,  9.0,
  //                             10.0, 10.0, 10.0, 11.0, 11.0, 11.0, 11.0, 12.0, 16.0, 2.0,  3.0,
  //                             5.0,  7.0,  9.0,  9.0,  10.0, 10.0, 10.0};
  //
  //  // testing dataset id:2
  //  std::vector<double> D = {0,  0,  0,  0,  0,   0,   0,   0,   0,   0,   25,  25,  25,
  //                           25, 25, 25, 25, 25,  25,  25,  50,  50,  50,  50,  50,  50,
  //                           50, 50, 50, 50, 100, 100, 100, 100, 100, 100, 100, 100, 100};
  //  std::vector<double> LS = {16, 9,  15, 14, 13, 9,  10, 14, 10, 11, 14, 9,  14,
  //                            9,  13, 12, 10, 10, 11, 14, 11, 11, 14, 11, 10, 11,
  //                            10, 15, 7,  14, 11, 14, 12, 13, 12, 14, 11, 8,  10};
  //  std::vector<double> I = {1, 1, 2, 3, 3, 0, 2, 2, 1, 2, 4, 5, 6, 2, 6, 3, 1, 2, 4, 3,
  //                           4, 5, 5, 4, 5, 4, 5, 6, 2, 4, 6, 6, 8, 7, 8, 6, 6, 5, 4};
  //  std::vector<double> LSC = {1, 1, 2, 3, 3, 0, 2, 2, 1, 2, 4, 5, 6, 2, 6, 3, 1, 2, 4, 3,
  //                             4, 5, 5, 4, 5, 4, 5, 6, 2, 4, 6, 6, 8, 7, 8, 6, 6, 5, 4};
  //
  // testing dataset id:3
  // std::vector<double> D = {0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
  //                          0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
  //                          0,    500,  500,  500,  500,  500,  500,  500,  500,  500,  500,  500,
  //                          500,  500,  500,  500,  500,  500,  500,  500,  500,  500,  500,  500,
  //                          500,  1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
  //                          1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
  //                          1000, 1000, 1000, 1000, 1000, 2000, 2000, 2000, 2000, 2000, 2000,
  //                          2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000,
  //                          2000, 2000, 2000, 2000, 2000, 2000, 2000, 5000, 5000, 5000, 5000,
  //                          5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000,
  //                          5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000};
  // std::vector<double> LS = {6, 6, 5, 6, 2, 7, 5, 6, 6, 6, 6, 5, 5, 1, 3, 3, 6, 2, 8, 2, 4,
  //                           7, 5, 7, 5, 6, 4, 2, 8, 4, 4, 7, 4, 6, 8, 3, 5, 7, 6, 2, 7, 7,
  //                           7, 5, 6, 8, 4, 8, 2, 1, 5, 6, 2, 7, 1, 3, 5, 2, 4, 8, 5, 8, 5,
  //                           8, 7, 4, 5, 7, 3, 3, 3, 7, 6, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6, 8,
  //                           7, 5, 5, 6, 5, 6, 6, 6, 6, 7, 7, 6, 9, 8, 6, 5, 4, 5, 6, 7, 6,
  //                           7, 1, 2, 6, 4, 3, 5, 7, 5, 4, 5, 6, 7, 4, 6, 5, 1, 6, 3, 4};
  // std::vector<double> I = {0, 2, 2, 0, 0, 3, 0, 2, 1, 0, 1, 1, 0, 0, 0, 1, 0, 2, 1, 1, 0,
  //                          0, 1, 1, 1, 1, 0, 0, 1, 2, 0, 3, 0, 2, 1, 0, 1, 0, 0, 0, 1, 2,
  //                          3, 1, 0, 1, 0, 0, 0, 0, 0, 3, 0, 0, 0, 1, 0, 1, 0, 2, 0, 3, 4,
  //                          0, 2, 1, 1, 2, 0, 1, 3, 1, 2, 0, 1, 1, 3, 2, 4, 2, 3, 4, 0, 5,
  //                          1, 1, 3, 2, 5, 0, 2, 5, 6, 3, 5, 5, 6, 7, 3, 1, 4, 4, 0, 6, 0,
  //                          7, 1, 2, 6, 4, 3, 4, 5, 5, 3, 5, 5, 7, 2, 4, 3, 1, 5, 1, 4};
  // std::vector<double> LSC =
  // {6.0,  13.0, 6.0,  7.0,  11.0, 12.0, 14.0, 12.0, 9.0,  9.0,  11.0, 7.0,
  //                            9.0,  13.0, 4.0,  13.0, 9.0,  6.0,  9.0,  5.0,  6.0,  7.0,  11.0, 9.0,
  //                            14.0, 7.0,  4.0,  2.0,  13.0, 12.0, 7.0,  7.0,  9.0,  9.0,  8.0,  14.0,
  //                            9.0,  7.0,  9.0,  9.0,  10.0, 7.0,  11.0, 5.0,  10.0, 10.0, 8.0,  14.0,
  //                            10.0, 13.0, 10.0, 14.0, 14.0, 8.0,  6.0,  3.0,  10.0, 9.0,  10.0, 9.0,
  //                            14.0, 8.0,  14.0, 13.0, 10.0, 14.0, 6.0,  12.0, 9.0,  7.0,  14.0, 11.0,
  //                            14.0, 9.0,  7.0,  11.0, 10.0, 14.0, 10.0, 8.0,  13.0, 11.0, 10.0, 13.0,
  //                            7.0,  6.0,  14.0, 6.0,  7.0,  14.0, 8.0,  12.0, 6.0,  10.0, 13.0, 12.0,
  //                            10.0, 11.0, 8.0,  6.0,  5.0,  14.0, 14.0, 13.0, 12.0, 9.0,  4.0,  13.0,
  //                            6.0,  12.0, 3.0,  7.0,  7.0,  10.0, 14.0, 11.0, 14.0, 7.0,  4.0,  8.0,
  //                            10.0, 13.0, 7.0,  10.0, 4.0};

  //  // testing dataset id:4
  //  std::vector<double> D = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
  //                           1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
  //                           3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};
  //  std::vector<double> LS = {17, 14, 9,  10, 12, 10, 10, 11, 10, 6,  7,  7,  15, 14, 10, 13, 12,
  //  11,
  //                            12, 14, 15, 12, 11, 12, 12, 14, 12, 11, 13, 13, 12, 11, 12, 14, 9,
  //                            12, 9,  17, 13, 10, 12, 10, 8,  10, 4,  9,  9,  6,  4,  7,  13, 10,
  //                            13, 11};
  //  std::vector<double> I = {3, 2,  2, 3, 2,  2,  2, 2,  2, 2, 2, 2, 2, 2, 2, 2, 2,  2,
  //                           2, 2,  3, 2, 3,  2,  4, 2,  4, 2, 3, 2, 2, 2, 2, 2, 2,  2,
  //                           2, 17, 9, 8, 12, 10, 8, 10, 4, 6, 9, 5, 4, 7, 3, 7, 10, 4};
  //  std::vector<double> LSC = {31.9, 29.2, 25.3, 26.3, 26.7, 27.8, 25.7, 27.5, 28.9, 27.4, 27.1,
  //                             26,   30.2, 28.9, 26.4, 24.9, 24.2, 26.9, 24.2, 26.9, 28.1, 30.6,
  //                             26.1, 29.3, 31.1, 28.1, 24.3, 25.6, 25.3, 27,   26.8, 26.7, 26.3,
  //                             26.8, 26.7, 26.8, 30.9, 30.1, 26.2, 25.6, 24.9, 25.5, 26.1, 25.7,
  //                             26.8, 28.1, 27.4, 25.2, 26.4, 26,   26.2, 27.5, 28.6, 25.5};
  //
  //  // testing dataset id:5;
  //  std::vector<double> D = {0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
  //                           0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
  //                           0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
  //                           0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
  //                           0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0, 1000,
  //                           1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
  //                           1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
  //                           1000, 1000, 1000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000,
  //                           2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000,
  //                           2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000,
  //                           2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000,
  //                           5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000,
  //                           5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000,
  //                           5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000,
  //                           5000, 5000, 5000, 5000, 5000, 5000, 5000};
  //  std::vector<double> LS = {6, 6, 5, 6, 2, 7, 5, 6, 6, 6, 6, 5, 5, 1, 3, 3, 6, 2, 8, 2, 4, 7, 5,
  //  7,
  //                            5, 2, 6, 7, 7, 8, 4, 7, 5, 6, 7, 6, 2, 5, 1, 7, 7, 6, 8, 7, 7, 4, 3,
  //                            6, 3, 4, 7, 5, 6, 6, 8, 5, 8, 6, 7, 1, 5, 6, 2, 7, 1, 3, 5, 2, 4, 8,
  //                            5, 8, 5, 8, 7, 4, 5, 7, 3, 3, 3, 7, 6, 7, 7, 7, 7, 6, 6, 6, 6, 6, 6,
  //                            8, 7, 5, 5, 6, 5, 6, 6, 6, 6, 7, 7, 6, 9, 8, 6, 9, 7, 8, 7, 7, 6, 7,
  //                            5, 2, 6, 7, 7, 4, 3, 4, 7, 6, 5, 4, 5, 6, 7, 6, 7, 1, 2, 6, 4, 3, 5,
  //                            7, 5, 4, 5, 6, 7, 4, 6, 5, 1, 6, 3, 4, 6, 2, 6, 6, 6, 6, 7, 3, 6, 4,
  //                            1, 8, 3, 4};
  //  std::vector<double> I = {0, 6, 3, 0, 0, 5, 0, 3, 1, 0, 3, 1, 0, 0, 0, 1, 0, 2, 3, 1, 0, 0, 1,
  //  1,
  //                           1, 0, 1, 3, 2, 1, 1, 1, 0, 3, 0, 2, 0, 1, 0, 5, 3, 5, 8, 2, 2, 0, 0,
  //                           3, 1, 1, 1, 0, 2, 5, 4, 1, 2, 3, 2, 0, 0, 3, 0, 0, 0, 1, 0, 2, 0, 2,
  //                           0, 8, 4, 0, 2, 3, 1, 2, 2, 1, 3, 2, 2, 2, 3, 3, 3, 2, 4, 2, 3, 4, 0,
  //                           5, 1, 1, 3, 2, 5, 0, 2, 5, 6, 3, 5, 5, 6, 7, 3, 0, 4, 0, 1, 4, 3, 1,
  //                           3, 2, 5, 5, 4, 1, 3, 2, 0, 3, 1, 4, 4, 0, 6, 0, 7, 1, 2, 6, 4, 3, 4,
  //                           5, 5, 3, 5, 5, 7, 2, 4, 3, 1, 5, 1, 4, 6, 1, 3, 5, 4, 2, 5, 1, 6, 4,
  //                           1, 0, 3, 4};
  //  std::vector<double> LSC = {6, 6, 5, 6, 2, 7, 5, 6, 6, 6, 6, 5, 5, 1, 3, 3, 6, 2, 8, 2, 4, 7,
  //  5, 7,
  //                             5, 2, 6, 7, 7, 8, 4, 7, 5, 6, 7, 6, 2, 5, 1, 7, 7, 6, 8, 7, 7, 4,
  //                             3, 6, 3, 4, 7, 5, 6, 6, 8, 5, 8, 6, 7, 1, 5, 6, 2, 7, 1, 3, 5, 2,
  //                             4, 8, 5, 8, 5, 8, 7, 4, 5, 7, 3, 3, 3, 7, 6, 7, 7, 7, 7, 6, 6, 6,
  //                             6, 6, 6, 8, 7, 5, 5, 6, 5, 6, 6, 6, 6, 7, 7, 6, 9, 8, 6, 9, 7, 8,
  //                             7, 7, 6, 7, 5, 2, 6, 7, 7, 4, 3, 4, 7, 6, 5, 4, 5, 6, 7, 6, 7, 1,
  //                             2, 6, 4, 3, 5, 7, 5, 4, 5, 6, 7, 4, 6, 5, 1, 6, 3, 4, 6, 2, 6, 6,
  //                             6, 6, 7, 3, 6, 4, 1, 8, 3, 4};
  //
  //  // testing dataset id:6;
  //  std::vector<double> D = {0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  //                           0,   0,   0,   0,   0,   0,   0,   0,   0,   100, 100, 100, 100, 100,
  //                           100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100,
  //                           100, 100, 100, 100, 100, 100, 250, 250, 250, 250, 250, 250, 250, 250,
  //                           250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250,
  //                           250};
  //  std::vector < double LS = {15, 17, 14, 14, 15, 15, 18, 12, 15, 15, 15, 19, 15, 16, 18, 18, 19,
  //  13,
  //                             16, 13, 8,  14, 18, 15, 14, 14, 13, 16, 16, 15, 14, 16, 14, 14, 15,
  //                             16, 12, 16, 16, 1,  14, 18, 16, 16, 15, 15, 15, 12, 18, 16, 16, 15,
  //                             15, 17, 14, 15, 13, 15, 17, 16, 16, 11, 15, 12, 6,  6,  2,  18, 18,
  //                             12};
  //  std::vector<double> I = {0, 0, 0, 0, 1, 0, 2, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1,
  //  0,
  //                           0, 2, 1, 1, 0, 0, 1, 0, 0, 3, 1, 1, 3, 2, 0, 0, 0, 0, 0, 1, 0, 1, 0,
  //                           0, 1, 2, 5, 1, 2, 0, 1, 0, 0, 0, 1, 1, 1, 3, 0, 1, 2, 0, 0, 2, 2, 4};
  //  std::vector<double> LSC = {15, 17, 14, 14, 15, 15, 18, 12, 15, 15, 15, 19, 15, 16, 18, 18, 19,
  //  13,
  //                             16, 13, 8,  14, 18, 15, 14, 14, 13, 16, 16, 15, 14, 16, 14, 14, 15,
  //                             16, 12, 16, 16, 1,  14, 18, 16, 16, 15, 15, 15, 12, 18, 16, 16, 15,
  //                             15, 17, 14, 15, 13, 15, 17, 16, 16, 11, 15, 12, 6,  6,  2,  18, 18,
  //                             12};
  //
  //  // testing dataset id:7;
  //  std::vector<double> D = {0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  //                           0,   0,   0,   0,   0,   0,   0,   0,   0,   100, 100, 100, 100, 100,
  //                           100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100,
  //                           100, 100, 100, 100, 100, 100, 250, 250, 250, 250, 250, 250, 250, 250,
  //                           250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250,
  //                           500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500};
  //  std::vector < double LS = {15, 17, 14, 14, 15, 15, 18, 12, 15, 15, 15, 19, 15, 16, 18, 18, 19,
  //                             13, 16, 13, 8,  14, 18, 15, 14, 14, 13, 16, 16, 15, 14, 16, 14, 14,
  //                             15, 16, 12, 16, 16, 1,  14, 18, 16, 16, 15, 15, 15, 12, 18, 16, 16,
  //                             15, 15, 17, 14, 15, 13, 15, 17, 16, 16, 11, 15, 12, 6,  6,  2,  18,
  //                             18, 12, 5,  12, 5,  15, 12, 16, 9,  6,  6,  11, 2};
  //  std::vector<double> I = {0, 0, 0, 0, 1, 0, 2, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0,
  //                           1, 1, 0, 0, 2, 1, 1, 0, 0, 1, 0, 0, 3, 1, 1, 3, 2, 0, 0, 0, 0,
  //                           0, 1, 0, 1, 0, 0, 1, 2, 5, 1, 2, 0, 1, 0, 0, 0, 1, 1, 1, 3, 0,
  //                           1, 2, 0, 0, 2, 2, 4, 0, 0, 1, 2, 1, 0, 1, 0, 1, 0, 0};
  //  std::vector<double> LSC = {170, 160, 147, 153, 158, 153, 168, 165, 164, 166, 149, 174, 156,
  //  160,
  //                             158, 161, 166, 172, 181, 177, 141, 144, 157, 161, 159, 153, 146,
  //                             167, 150, 159, 152, 165, 166, 158, 168, 143, 148, 177, 154, 153,
  //                             179, 171, 180, 170, 165, 157, 164, 162, 159, 160, 151, 141, 179,
  //                             150, 153, 175, 146, 161, 167, 165, 166, 162, 157, 153, 158, 166,
  //                             167, 146, 164, 155, 161, 158, 181, 159, 151, 152, 166, 176, 165,
  //                             144, 144};
  //
  //  // testing dataset id:8;
  //  std::vector<double> D = {0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  //                           125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 250, 250,
  //                           250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250,
  //                           500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 625,
  //                           625, 625, 625, 625, 625, 625, 625, 625, 625};
  //
  //  std::vector < double LS = {4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  //                             4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4,
  //                             4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4};
  //  std::vector<double> I = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 1, 0, 0,
  //                           1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 3, 3, 0, 0, 1, 4,
  //                           2, 1, 2, 2, 4, 1, 3, 2, 2, 4, 3, 3, 4, 2, 4, 4, 3, 4, 4};
  //  std::vector<double> LSC = {12, 17, 12, 15, 18, 13, 16, 11, 13, 14, 11, 16, 12, 14, 14,
  //                             12, 15, 15, 13, 16, 16, 18, 11, 12, 14, 14, 16, 14, 14, 16,
  //                             10, 12, 15, 14, 13, 10, 14, 17, 15, 17, 15, 16, 11, 14, 15,
  //                             11, 13, 11, 12, 16, 15, 11, 11, 15, 12, 13, 19};
  //
  //  // testing dataset id:9;
  //  std::vector<double> D = {0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
  //                           125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 250, 250,
  //                           250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250,
  //                           500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 625,
  //                           625, 625, 625, 625, 625, 625, 625, 625, 625};
  //  std::vector < double LS = {4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  //                             4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4,
  //                             4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4};
  //  std::vector<double> I = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 1, 0, 0,
  //                           1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 3, 3, 0, 0, 1, 4,
  //                           2, 1, 2, 2, 4, 1, 3, 2, 2, 4, 3, 3, 4, 2, 4, 4, 3, 4, 4};
  //  std::vector<double> LSC = {4, 5, 3,  6, 8,  8, 5, 7,  10, 6, 5, 12, 6, 7, 9, 7,  6, 10, 7,
  //                             8, 9, 11, 5, 10, 4, 9, 12, 6,  8, 9, 3,  7, 7, 9, 7,  4, 7,  12,
  //                             6, 9, 7,  7, 5,  6, 5, 6,  11, 7, 4, 9,  7, 8, 8, 10, 8, 9,  8};
  //
  //  // testing dataset id:10
  //  std::vector<double> D = {0,  0,  0,  0,  0,   0,   0,   0,   0,   0,   25,  25,  25,
  //                           25, 25, 25, 25, 25,  25,  25,  50,  50,  50,  50,  50,  50,
  //                           50, 50, 50, 50, 100, 100, 100, 100, 100, 100, 100, 100, 100};
  //  std::vector<double> LS = {16, 9,  15, 14, 13, 9,  10, 14, 10, 11, 14, 9,  14,
  //                            9,  13, 12, 10, 10, 11, 14, 11, 11, 14, 11, 10, 11,
  //                            10, 15, 7,  14, 11, 14, 12, 13, 12, 14, 11, 8,  10};
  //  std::vector<double> I = {1, 1, 2, 3, 3, 0, 2, 2, 1, 2, 4, 5, 6, 2, 6, 3, 1, 2, 4, 3,
  //                           4, 5, 5, 4, 5, 4, 5, 6, 2, 4, 6, 6, 8, 7, 8, 6, 6, 5, 4};
  //  std::vector<double> LSC = {16.0, 9.0,  15.0, 14.0, 13.0, 9.0,  10.0, 14.0, 10.0, 11.0,
  //                             14.0, 9.0,  14.0, 9.0,  13.0, 12.0, 10.0, 10.0, 11.0, 14.0,
  //                             11.0, 11.0, 14.0, 11.0, 10.0, 11.0, 10.0, 15.0, 7.0,  14.0,
  //                             11.0, 14.0, 12.0, 13.0, 12.0, 14.0, 11.0, 8.0,  10.0};

  // parameter counting test
  std::vector<double> D = {1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4};
  std::vector<double> LS = {5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0,
                            5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0};
  std::vector<double> I = {1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4};
  std::vector<double> LSC = {1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4};

  pyAnal.doses = D;
  pyAnal.litterSize = LS;
  pyAnal.incidence = I;
  pyAnal.lsc = LSC;

  pyAnal.LSC_type = 1;  // 1 = Overall Mean; 2 = control group mean; 0 = do not use LSC
  pyAnal.ILC_type = 1;  // 1 = estimate intralitter; assume 0 otherwise
  pyAnal.BMD_type = 1;  // 1 = extra; added otherwise
  pyAnal.estBackground = true;
  pyAnal.BMR = 0.1;
  pyAnal.alpha = 0.05;
  pyAnal.numBootRuns = 3;
  pyAnal.iterations = 1000;
  pyAnal.seed = BMDS_MISSING;
  pyAnal.prior_cols = 2;

  int numDoseGroups = 4;

  int Nobs = pyAnal.doses.size();
  pyAnal.prior.resize(pyAnal.prior_cols * (5 + numDoseGroups), 0.0);
  if (pyAnal.model == nlogistic) {
    pyAnal.prior = getNlogisticPrior(numDoseGroups, pyAnal.prior_cols, isRestricted);
  } else if (pyAnal.model == nctr) {
    pyAnal.prior = getNctrPrior(numDoseGroups, pyAnal.prior_cols, isRestricted);
  } else {
    std::cout << "Incorrect model selected for testing" << std::endl;
  }

  struct python_nested_result pyRes;
  pyAnal.parms = 5 + numDoseGroups;
  pyRes.model = pyAnal.model;
  pyRes.nparms = pyAnal.parms;

  struct BMDS_results bmdsRes;

  // set all parms as unbounded initially
  for (int i = 0; i < pyAnal.parms; i++) {
    bmdsRes.bounded.push_back(false);
    bmdsRes.stdErr.push_back(BMDS_MISSING);
    bmdsRes.lowerConf.push_back(BMDS_MISSING);
    bmdsRes.upperConf.push_back(BMDS_MISSING);
  }
  bmdsRes.BMD = -9999.0;
  bmdsRes.BMDU = -9999.0;
  bmdsRes.BMDL = -9999.0;
  bmdsRes.AIC = -9999.0;
  pyRes.bmdsRes = bmdsRes;

  struct nestedLitterData litter;
  std::vector<double> litterDose(Nobs);
  std::vector<double> litterLSC(Nobs);
  std::vector<double> litterEP(Nobs);
  std::vector<double> litterSize(Nobs);
  std::vector<double> litterEx(Nobs);
  std::vector<int> litterObs(Nobs);
  std::vector<double> litterSR(Nobs);
  litter.dose = litterDose;
  litter.LSC = litterLSC;
  litter.estProb = litterEP;
  litter.litterSize = litterSize;
  litter.expected = litterEx;
  litter.observed = litterObs;
  litter.SR = litterSR;
  pyRes.litter = litter;

  struct nestedReducedData redData;
  std::vector<double> redDose(numDoseGroups);  // size = numRows
  std::vector<double> redPA(numDoseGroups);    // estimate of proportion affected
  std::vector<double> redLC(numDoseGroups);    // reduced data lower confidence limit
  std::vector<double> redUC(numDoseGroups);    // reduced data upper confidence limit
  redData.dose = redDose;                      // size = numRows
  redData.propAffect = redPA;                  // estimate of proportion affected
  redData.lowerConf = redLC;
  redData.upperConf = redUC;
  pyRes.reduced = redData;

  struct nestedSRData srData;
  pyRes.srData = srData;

  struct nestedBootstrap bootData;
  int numBootRuns = 3;
  std::vector<double> pVal(numBootRuns + 1);
  std::vector<double> perc50(numBootRuns + 1);
  std::vector<double> perc90(numBootRuns + 1);
  std::vector<double> perc95(numBootRuns + 1);
  std::vector<double> perc99(numBootRuns + 1);
  bootData.pVal = pVal;
  bootData.perc50 = perc50;
  bootData.perc90 = perc90;
  bootData.perc95 = perc95;
  bootData.perc99 = perc99;
  pyRes.boot = bootData;

  pythonBMDSNested(&pyAnal, &pyRes);

  std::cout << "INPUT" << std::endl;
  printBmdsStruct(&pyAnal);
  std::cout << "OUTPUT" << std::endl;
  printBmdsStruct(&pyRes);
  // printNestedModResult(&pyAnal, &pyRes, showResultsOverride);
}

void printNestedModResult(
    struct python_nested_analysis *pyAnal, struct python_nested_result *pyRes,
    bool showResultsOverride
) {
  printf("tlink pyRes.validResult = %s\n", pyRes->bmdsRes.validResult ? "valid" : "invalid");
  if (pyRes->bmdsRes.validResult || showResultsOverride) {
    // std::cout<<"Valid Result"<<std::endl;
    printf("\nBenchmark Dose\n");
    printf("BMD: %f\n", pyRes->bmdsRes.BMD);
    printf("BMDL: %f\n", pyRes->bmdsRes.BMDL);
    printf("BMDU: %f\n", pyRes->bmdsRes.BMDU);
    printf("LL: %f\n", pyRes->LL);
    printf("AIC: %f\n", pyRes->bmdsRes.AIC);
    printf("P-value: %f\n", pyRes->combPVal);
    printf("Chi^2: %f\n", pyRes->bmdsRes.chisq);

    printf("\nModel Parameters\n");
    printf("# of parms: %d\n", pyRes->nparms);
    printf("DOF: %f\n", pyRes->model_df);
    //      printf("parm, estimate, bounded, std.err., lower conf, upper conf\n");
    for (int i = 0; i < pyRes->nparms; i++) {
      printf(
          "%d, %.10f, %s, %f\n", i, pyRes->parms[i], pyRes->bmdsRes.bounded[i] ? "true" : "false",
          pyRes->bmdsRes.stdErr[i]
      );
    }
  }
  //   printf("\ncov matrix\n");
  //   for (int i=0; i<pyAnal->parms*pyAnal->parms; i++){
  //     printf("%d, %f\n", i, pyRes->cov[i]);
  //   }

  std::cout << "----SROI Data----" << std::endl;
  std::cout << "min SR:" << pyRes->srData.minSR << std::endl;
  std::cout << "min |SR|:" << pyRes->srData.minAbsSR << std::endl;
  std::cout << "avg SR:" << pyRes->srData.avgSR << std::endl;
  std::cout << "avg |SR|:" << pyRes->srData.avgAbsSR << std::endl;
  std::cout << "max SR:" << pyRes->srData.maxSR << std::endl;
  std::cout << "max |SR|:" << pyRes->srData.maxAbsSR << std::endl;

  std::cout << "----Litter Data----" << std::endl;
  std::cout << "Dose\tLSC\tEstProb\t\tLS\tExp\tObs\tSR" << std::endl;
  for (int i = 0; i < pyRes->litter.dose.size(); i++) {
    std::cout << pyRes->litter.dose[i] << "\t" << pyRes->litter.LSC[i] << "\t"
              << pyRes->litter.estProb[i] << "\t" << pyRes->litter.litterSize[i] << "\t"
              << pyRes->litter.expected[i] << "\t" << pyRes->litter.observed[i] << "\t"
              << pyRes->litter.SR[i] << std::endl;
  }
  std::cout << "chiSq:" << pyRes->litter.chiSq << std::endl;

std:
  cout << "----Bootstrap Data----" << std::endl;
  int numRuns = pyRes->boot.pVal.size();
  for (int i = 0; i < numRuns; i++) {
    std::cout << "i:" << i << ", pval:" << pyRes->boot.pVal[i] << ", 50th:" << pyRes->boot.perc50[i]
              << ", 90th:" << pyRes->boot.perc90[i] << ", 95th:" << pyRes->boot.perc95[i]
              << ", 99th:" << pyRes->boot.perc99[i] << std::endl;
  }
}

void Nlogist_probs_test() {
  std::vector<double> Xi = {0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
                            0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
                            0,    500,  500,  500,  500,  500,  500,  500,  500,  500,  500,  500,
                            500,  500,  500,  500,  500,  500,  500,  500,  500,  500,  500,  500,
                            500,  1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
                            1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
                            1000, 1000, 1000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000,
                            2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000,
                            2000, 2000, 2000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000,
                            5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000,
                            5000, 5000, 5000, 5000, 5000};
  std::vector<double> Ls = {4,  5,  6,  6,  6,  6,  7,  7,  7,  9,  9,  9,  9,  9,  9,  11, 11, 11,
                            12, 12, 13, 13, 13, 14, 14, 2,  4,  5,  7,  7,  7,  7,  7,  8,  8,  9,
                            9,  9,  9,  9,  10, 10, 10, 10, 11, 12, 13, 14, 14, 3,  6,  6,  7,  7,
                            8,  8,  9,  9,  9,  9,  10, 10, 10, 10, 11, 12, 13, 13, 14, 14, 14, 14,
                            14, 14, 14, 6,  6,  6,  7,  7,  8,  8,  8,  10, 10, 10, 10, 10, 11, 11,
                            11, 12, 12, 13, 13, 13, 14, 14, 14, 3,  4,  4,  4,  5,  6,  6,  7,  7,
                            7,  7,  8,  9,  10, 10, 10, 11, 12, 12, 13, 13, 13, 14, 14, 14, 14};

  std::vector<double> p = {0.054599, -10.777510, 0.165150, -0.059984, 1.437543,
                           0.015216, 0.047678,   0.247516, 0.254625,  0.551903};
  // only care about Spec[0] and Spec[2]
  std::vector<bool> Spec = {false, false, false, false, false, false, false, false, false, false};
  bool compgrad = false;

  int Nobs = Xi.size();
  std::vector<double> probs(Nobs);
  std::vector<std::vector<double>> gradij(Nobs, std::vector<double>(5));
  struct nestedObjData objData;
  objData.isBMDL = false;
  objData.smax = 14.0;
  objData.smin = 2.0;
  objData.Ls = Ls;
  objData.Xi = Xi;
  objData.sijfixed = 9.28;
  objData.riskType = 1;
  objData.BMR = 0.1;
  objData.tD = 575.903266;
  objData.Spec = Spec;

  Nlogist_probs(probs, p, compgrad, gradij, &objData);
  for (int i = 0; i < Nobs; i++) {
    std::cout << "i:" << i << ", probs:" << probs[i] << std::endl;
  }
}

void Nlogist_lk_test() {
  std::vector<double> Xi = {0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
                            0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
                            0,    500,  500,  500,  500,  500,  500,  500,  500,  500,  500,  500,
                            500,  500,  500,  500,  500,  500,  500,  500,  500,  500,  500,  500,
                            500,  1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
                            1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000,
                            1000, 1000, 1000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000,
                            2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000, 2000,
                            2000, 2000, 2000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000,
                            5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000,
                            5000, 5000, 5000, 5000, 5000};
  std::vector<double> Ls = {4,  5,  6,  6,  6,  6,  7,  7,  7,  9,  9,  9,  9,  9,  9,  11, 11, 11,
                            12, 12, 13, 13, 13, 14, 14, 2,  4,  5,  7,  7,  7,  7,  7,  8,  8,  9,
                            9,  9,  9,  9,  10, 10, 10, 10, 11, 12, 13, 14, 14, 3,  6,  6,  7,  7,
                            8,  8,  9,  9,  9,  9,  10, 10, 10, 10, 11, 12, 13, 13, 14, 14, 14, 14,
                            14, 14, 14, 6,  6,  6,  7,  7,  8,  8,  8,  10, 10, 10, 10, 10, 11, 11,
                            11, 12, 12, 13, 13, 13, 14, 14, 14, 3,  4,  4,  4,  5,  6,  6,  7,  7,
                            7,  7,  8,  9,  10, 10, 10, 11, 12, 12, 13, 13, 13, 14, 14, 14, 14};

  std::vector<double> Yp = {0, 1, 0, 0, 2, 2, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 3, 2, 0,
                            1, 2, 0, 1, 0, 0, 1, 1, 0, 2, 0, 3, 0, 1, 0, 0, 2, 1, 0, 0, 1,
                            1, 0, 3, 2, 1, 0, 0, 1, 0, 1, 1, 1, 3, 0, 2, 1, 0, 0, 0, 0, 0,
                            2, 1, 2, 0, 0, 4, 0, 0, 3, 2, 3, 1, 6, 2, 1, 5, 1, 3, 2, 2, 6,
                            3, 0, 4, 3, 7, 4, 1, 5, 5, 5, 3, 5, 0, 3, 2, 3, 4, 1, 2, 4, 6,
                            1, 7, 5, 4, 5, 4, 7, 1, 3, 5, 5, 4, 0, 1, 2, 6, 5, 3, 0, 4};
  std::vector<double> Yn = {3, 1, 6, 4, 3, 0, 7, 6, 4, 5, 6, 6, 6, 5, 7, 4, 5, 2, 4, 4, 1,
                            2, 4, 5, 4, 2, 4, 4, 5, 4, 5, 7, 4, 4, 7, 2, 6, 4, 4, 4, 6, 6,
                            7, 2, 4, 2, 7, 8, 3, 2, 1, 4, 2, 6, 5, 7, 6, 1, 7, 3, 4, 5, 5,
                            5, 6, 5, 1, 8, 1, 2, 5, 3, 4, 0, 3, 0, 4, 4, 0, 6, 3, 4, 4, 3,
                            4, 6, 2, 4, 1, 2, 6, 1, 1, 3, 3, 2, 6, 2, 4, 0, 0, 0, 2, 0, 0,
                            4, 0, 2, 1, 1, 2, 0, 2, 2, 0, 0, 0, 6, 0, 0, 1, 1, 1, 6, 1};
  std::vector<int> Xg = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                         2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3,
                         3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                         4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5,
                         5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};

  std::vector<double> p = {0.054599, -10.777510, 0.165150, -0.059984, 1.437543,
                           0.015216, 0.047678,   0.247516, 0.254625,  0.551903};
  // only care about Spec[0] and Spec[2]
  std::vector<bool> Spec = {false, false, false, false, false, false, false, false, false, false};
  bool compgrad = false;

  int Nobs = Xi.size();
  std::vector<double> probs(Nobs);
  std::vector<std::vector<double>> gradij(Nobs, std::vector<double>(5));
  struct nestedObjData objData;
  objData.isBMDL = false;
  objData.smax = 14.0;
  objData.smin = 2.0;
  objData.Ls = Ls;
  objData.Xi = Xi;
  objData.Yp = Yp;
  objData.Yn = Yn;
  objData.Xg = Xg;
  objData.sijfixed = 9.28;
  objData.riskType = 1;
  objData.BMR = 0.1;
  objData.tD = 575.903266;
  objData.Spec = Spec;

  Nlogist_probs(probs, p, compgrad, gradij, &objData);
  for (int i = 0; i < Nobs; i++) {
    std::cout << "i:" << i << ", probs:" << probs[i] << std::endl;
  }
  double lk = Nlogist_lk(p, &objData);
  std::cout << "actual lk:" << lk << std::endl;
}

// void Nlogist_lk_test(){
//
//         std::vector<double> Xi =
//         {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,2000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000,5000};
//         std::vector<double> Ls =
//         {4,5,6,6,6,6,7,7,7,9,9,9,9,9,9,11,11,11,12,12,13,13,13,14,14,2,4,5,7,7,7,7,7,8,8,9,9,9,9,9,10,10,10,10,11,12,13,14,14,3,6,6,7,7,8,8,9,9,9,9,10,10,10,10,11,12,13,13,14,14,14,14,14,14,14,6,6,6,7,7,8,8,8,10,10,10,10,10,11,11,11,12,12,13,13,13,14,14,14,3,4,4,4,5,6,6,7,7,7,7,8,9,10,10,10,11,12,12,13,13,13,14,14,14,14};
//
//         std::vector<double> Yp =
//         {0,1,0,0,2,2,0,0,1,1,1,0,0,0,1,1,1,0,3,2,0,1,2,0,1,0,0,1,1,0,2,0,3,0,1,0,0,2,1,0,0,1,1,0,3,2,1,0,0,1,0,1,1,1,3,0,2,1,0,0,0,0,0,2,1,2,0,0,4,0,0,3,2,3,1,6,2,1,5,1,3,2,2,6,3,0,4,3,7,4,1,5,5,5,3,5,0,3,2,3,4,1,2,4,6,1,7,5,4,5,4,7,1,3,5,5,4,0,1,2,6,5,3,0,4};
//	std::vector<double> Yn =
//{3,1,6,4,3,0,7,6,4,5,6,6,6,5,7,4,5,2,4,4,1,2,4,5,4,2,4,4,5,4,5,7,4,4,7,2,6,4,4,4,6,6,7,2,4,2,7,8,3,2,1,4,2,6,5,7,6,1,7,3,4,5,5,5,6,5,1,8,1,2,5,3,4,0,3,0,4,4,0,6,3,4,4,3,4,6,2,4,1,2,6,1,1,3,3,2,6,2,4,0,0,0,2,0,0,4,0,2,1,1,2,0,2,2,0,0,0,6,0,0,1,1,1,6,1};
//	std::vector<int> Xg =
//{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5};
//
//         std::vector<double> p =
//         {0.054599,-14.099836,0.165150,-0.059984,1.437543,0.015216,0.047678,0.247516,0.254625,0.551903};
//         //only care about Spec[0] and Spec[2]
//         std::vector<bool> Spec = {false, false, false, false, false, false, false, false, false,
//         false}; bool compgrad = false;
//
//         int Nobs = Xi.size();
//         std::vector<double> probs(Nobs);
//         std::vector<std::vector<double>> gradij(Nobs, std::vector<double> (5));
//         struct nestedObjData objData;
//         objData.isBMDL = false;
//         objData.smax = 14.0;
//         objData.smin = 2.0;
//         objData.Ls = Ls;
//         objData.Xi = Xi;
//	objData.Yp = Yp;
//	objData.Yn = Yn;
//	objData.Xg = Xg;
//         objData.Spec = Spec;
//         objData.sijfixed = 9.28;
//         objData.riskType = 1;
//         objData.BMR = 0.1;
//         objData.tD = 575.903266;
//
//
//         double exp_lk = -337.540036;
//         double lk = Nlogist_lk(p, &objData);
//
//         std::cout<<"expected lk:"<<exp_lk<<std::endl;
//	std::cout<<"actual lk:"<<lk<<std::endl;
//
// }

void runTestMTInequalityConstraint() {
  // USER INPUT
  // double target = 565.17482044575650;
  // double target = -104.341181352210455;
  double target = -103.61387669479440;

  const std::vector<double> x = {-2.8269246586783039, 1.3081031558978955E-021,
                                 0.69271407326031031, 1.1269537527262939E-021,
                                 -0.0000000000000000, 0.54665428796084770,
                                 0.16803414281323573, -0.0000000000000000,
                                 0.40754846974911302, 2.0785118843572268};
  // const std::vector<double> x =
  // {-2.8269246587415990,6.5404535862508677E-022,0.69271208814523166,5.6347646693637056E-022,-0.0000000000000000,0.54665480790790499,0.16803682933352726,1.9124921929208785E-021,0.40755047235337671,2.0785001207417775};

  std::vector<int> degree = {2, 2, 2};
  // data
  std::vector<double> doses1 = {0, 50, 100, 200, 400};
  std::vector<double> Y1 = {0, 1, 2, 10, 19};
  std::vector<double> n_group1 = {20, 20, 20, 20, 20};
  std::vector<double> doses2 = {0, 50, 100, 200, 400};
  std::vector<double> Y2 = {0, 1, 2, 4, 11};
  std::vector<double> n_group2 = {20, 20, 20, 20, 20};
  std::vector<double> doses3 = {0, 50, 100, 200, 400};
  std::vector<double> Y3 = {0, 2, 2, 6, 9};
  std::vector<double> n_group3 = {20, 20, 20, 20, 20};

  std::vector<std::vector<double>> doses;
  std::vector<std::vector<double>> Y;
  std::vector<std::vector<double>> n_group;
  doses.push_back(doses1);
  doses.push_back(doses2);
  doses.push_back(doses3);
  Y.push_back(Y1);
  Y.push_back(Y2);
  Y.push_back(Y3);
  n_group.push_back(n_group1);
  n_group.push_back(n_group2);
  n_group.push_back(n_group3);
  int nT = doses.size();

  // END USER INPUT

  double maxDose = 0;
  for (int i = 0; i < nT; i++) {
    double tmpMax = *std::max_element(doses[i].begin(), doses[i].end());
    if (tmpMax > maxDose) maxDose = tmpMax;
  }
  std::vector<double> grad(x.size());
  std::vector<int> nObs;
  struct msComboInEq ineq1;
  ineq1.nT = nT;
  ineq1.target = target;

  for (int i = 0; i < nT; i++) {
    std::vector<double> scaledDose = doses[i];
    nObs.push_back(scaledDose.size());
    std::cout << "dataset " << i << std::endl;
    for (int j = 0; j < scaledDose.size(); j++) {
      scaledDose[j] /= maxDose;
      std::cout << "j:" << j << ", dose:" << scaledDose[j] << std::endl;
    }
    ineq1.doses.push_back(scaledDose);
    ineq1.Y.push_back(Y[i]);
    ineq1.n_group.push_back(n_group[i]);
  }
  ineq1.nObs = nObs;
  ineq1.degree = degree;
  for (int i = 0; i < nT; i++) {
    std::cout << "i:" << i << ", nObs:" << nObs[i] << std::endl;
  }
  std::cout << "target:" << target << std::endl;

  std::cout << std::setprecision(15);
  double ineqVal = myInequalityConstraint1(x, grad, &ineq1);
  std::cout << "ineqVal = " << ineqVal << std::endl;

  std::cout << "ineq grad:" << std::endl;
  for (int i = 0; i < grad.size(); i++) {
    std::cout << "i:" << i << ", grad:" << grad[i] << std::endl;
  }
}

void runTestMTEqualityConstraint() {
  //  //USER INPUT

  const std::vector<double> x = {-2.8269246586783039, 1.3081031558978955E-021,
                                 0.69271407326031031, 1.1269537527262939E-021,
                                 -0.0000000000000000, 0.54665428796084770,
                                 0.16803414281323573, -0.0000000000000000,
                                 0.40754846974911302, 2.0785118843572268};
  // const std::vector<double> x =
  // {-2.8269246587415990,6.5404535862508677E-022,0.69271208814523166,5.6347646693637056E-022,-0.0000000000000000,0.54665480790790499,0.16803682933352726,1.9124921929208785E-021,0.40755047235337671,2.0785001207417775};

  std::vector<int> degree = {2, 2, 2};
  int nT = degree.size();
  double bmr = 0.1;
  std::vector<double> grad(x.size());

  // END USER INPUT

  struct msComboEq eq1;
  eq1.bmr = bmr;
  eq1.nT = nT;
  eq1.degree = degree;

  std::cout << std::setprecision(18);
  double ineqVal = myEqualityConstraint(x, grad, &eq1);
  std::cout << "eqVal = " << ineqVal << std::endl;

  std::cout << "eq grad:" << std::endl;
  for (int i = 0; i < grad.size(); i++) {
    std::cout << "i:" << i << ", grad:" << grad[i] << std::endl;
  }
}
