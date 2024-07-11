//integration_tests.cpp
//
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "priors.h"
#include <string>
#include <iostream>
#include "integration_tests.h"

bool showResultsOverride = true;


int run_all_integrationTests(){
	runPythonDichoAnalysis();
	runPythonContAnalysis();
	runPythonMultitumorAnalysis();
	return 0;
}

void runPythonDichoAnalysis(){

  printf("Running dichotomous analysis\n");

///////////////////////////////
//USER INPUT
///////////////////////////////

  enum dich_model model = d_multistage;  //d_hill =1, d_gamma=2,d_logistic=3, d_loglogistic=4,
                                   //d_logprobit=5, d_multistage=6,d_probit=7,
                                   //d_qlinear=8,d_weibull=9
  int modelType = 1;       //1 = frequentist, 2 = bayesian
  bool restricted = true;  //only used for frequentist models
  int BMD_type = 1;        // 1 = extra ; added otherwise
  int degree = 3;  //for multistage only
  double BMR = 0.1;
  double alpha = 0.05;
///////////////////////////////
//dicho data - dose, N, incidence
///////////////////////////////
  //Dichotomous.dax Effect 1
  double D[] = {0,50, 100, 150, 200};
  double Y[] = {0, 5, 30, 65, 90};
  double N[] = {100, 100, 100, 100, 100};

  /////////////////////////////////////////////////
////END USER INPUT
////////////////////////////////////////////////////

  //struct dichotomous_analysis anal;
  struct python_dichotomous_analysis anal;

  int numDataRows = sizeof(D)/sizeof(D[0]);

  //check data array sizes for consistency
  size_t numElementsY = sizeof(Y)/sizeof(Y[0]);
  size_t numElementsN = sizeof(N)/sizeof(N[0]);
  if (numDataRows != numElementsY || numElementsY != numElementsN) {
    printf("Number of data elements are not consistent\nExiting Code\n");
    exit(-1);
  }


  //priors defined columnwise
  int prCols = 5;

  //define priors/parameter constraints
  int numParms;
  printf("model = %d\n",model);
  switch(model) {
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
       //numParms = 2 + degree;
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
    default :
      printf("error in numParms\n");
      return;

  }

  double *prior;
  if (modelType == 1) {
    //frequentist
    if (restricted) {
      switch(model) {
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
          if (degree == 1){
            prior = prRFreqMulti1;
          } else if (degree == 2){
            prior = prRFreqMulti2;
          } else if (degree == 3){
            prior = prRFreqMulti3;
          } else if (degree == 4){
            prior = prRFreqMulti4;
          } else if (degree == 5){
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
      //unrestricted
      switch(model) {
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
          if (degree == 1){
            prior = prUFreqMulti1;
          } else if (degree == 2){
            prior = prUFreqMulti2;
          } else if (degree == 3){
            prior = prUFreqMulti3;
          } else if (degree == 4){
            prior = prUFreqMulti4;
          } else if (degree == 5){
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
    //bayesian
    switch(model) {
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
          if (degree == 1){
            prior = prBayesianMulti1;
          } else if (degree == 2){
            prior = prBayesianMulti2;
          } else if (degree == 3){
            prior = prBayesianMulti3;
          } else if (degree == 4){
            prior = prBayesianMulti4;
          } else if (degree == 5){
            prior = prBayesianMulti5;
          }
          break;
        case d_probit:
          anal.model = d_probit;
          prior = prBayesianProbit;
          break;
        case d_qlinear:
          anal.model =d_qlinear;
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

  //declare analysis
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
  anal.prior.assign(prior, prior + anal.prior_cols*anal.parms);

  struct python_dichotomous_model_result res;
  res.model = anal.model;
  res.dist_numE = 200;
  res.nparms = anal.parms;


  pythonBMDSDicho(&anal, &res);

//  printDichoModResult(&anal, &res, showResultsOverride);

}




void runPythonContAnalysis(){

  printf("Running python continuous analysis\n");

  //char * bmdsVersion[32];

  //version(bmdsVersion);
  string bmdsVersion = version();
  //printf("Version: %s\n", bmdsVersion);
  std::cout << "Version: " << bmdsVersion << std::endl;

  bool isIncreasing;

///////////////////////////////
//USER INPUT
//////////////////////////////

  enum cont_model model = exp_5; //hill, exp_3, exp_5, power, funl, polynomial
  int modelType = 1;   //1 = frequentist, 2 = bayesian
  bool restricted = true;   //only used for frequentist models
  enum distribution dist = normal;  //normal, normal_ncv, log_normal
  bool detectAdvDir = true;  //if false then need to set isIncreasing
  //isIncreasing = true;

  int degree = 2; //for polynomial only

  double alpha = 0.05;
  double BMRF = 1.0; //1.0;
  int BMD_type = 2;  //1=absdev, 2 = stddev, 3 = reldev, 4 = pt, 5 = extra, 6 = hybrid_extra, 7 = hybrid_added   from src/include/cmodeldefs.h
////////////////////////////////////////////////
//cont data - suff stat: dose, Y, N, SD
//cont data - individual: dose, response
/////////////////////////////////////////////
  bool suffStat = true;

  //continuous1.dax
  double D[] = {0,25,50, 100, 200};
  double Y[] = {6.0, 5.2, 2.4, 1.1, 0.75};
  double N[] = {20, 20, 19, 20, 20};
  double SD[] = {1.2, 1.1, 0.81, 0.74, 0.66};
  //isIncreasing = false;
/////////////////////////////////////////////////
//END USER INPUT
///////////////////////////////////////////////////



  //struct continuous_analysis anal;
  struct python_continuous_analysis anal;
  int numDataRows = sizeof(D)/sizeof(D[0]);

  if (!detectAdvDir){
    anal.isIncreasing = isIncreasing;
  }

  //check data array sizes for consistency
  size_t numElementsY = sizeof(Y)/sizeof(Y[0]);
  if (suffStat){
    size_t numElementsN = sizeof(N)/sizeof(N[0]);
    size_t numElementsSD = sizeof(SD)/sizeof(SD[0]);
    if (numDataRows != numElementsY || numElementsY != numElementsN || numElementsN != numElementsSD) {
      printf("Number of data elements are not consistent\nExiting Code\n");
      exit(-1);
    }
  } else {
    if (numDataRows != numElementsY) {
      printf("Number of data elements are not consistent\nExiting Code\n");
      exit(-1);
    }
  }

  //priors defined columnwise
  int prCols = 5;


  //define priors/parameter constraints
  int numParms;
  printf("model = %d\n",model);
  switch(model) {
    case hill:
       numParms = 6;
       break;
    case exp_3:
       //numParms = 5;
       //break;
    case exp_5:
       numParms = 6;
       break;
    case power:
       numParms = 5;
       break;
    case funl:
       numParms = 0;  //FIX THIS
       break;
    case polynomial:
       numParms = 3 + degree;
       break;
    default :
      printf("error in numParms\n");
      return;

  }
  if (dist == normal || dist == log_normal){
    numParms -= 1;
  }

  printf("numParms = %d\n", numParms);
  //double* pr;
  double *prior;


  printf("starting priors\n");

  if (modelType == 1) {
    //frequentist
    if (restricted) {
      printf("choosing frequentist restricted priors\n");
      switch(model) {
        case hill:
          anal.model = hill;
          if (dist == normal || dist == log_normal) {
            //normal
            prior = prRFreqHillNormal;
          } else {
          //} else if (dist == normal_ncv){
            //normal NCV
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
          if (dist == normal || dist == log_normal){
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
          if (detectAdvDir){

            if(dist == normal || dist == log_normal){
              printf("using advDir auto normal or log_normal dist priors\n");
              if (degree == 1) {
                prior = prRFreqPoly1;
              } else if (degree == 2){
                prior = prRFreqPoly2;
              } else if (degree == 3){
                prior = prRFreqPoly3;
              } else if (degree == 4){
                prior = prRFreqPoly4;
              } else if (degree == 5){
                prior = prRFreqPoly5;
              } else{
                printf("poly restricted normal/lognormal degree error\n");
                return;
              }
            } else {
              printf("using advDir auto normal_ncv dist priors\n");
              if (degree == 1) {
                prior = prRFreqPoly1NCV;
              } else if (degree == 2){
                prior = prRFreqPoly2NCV;
              } else if (degree == 3){
                prior = prRFreqPoly3NCV;
              } else if (degree == 4){
                prior = prRFreqPoly4NCV;
              } else if (degree == 5){
                prior = prRFreqPoly5NCV;
              } else{
                printf("poly restricted normal NCV degree error\n");
                return;
              }
            }
          } else {

            if (anal.isIncreasing) {
              if(dist == normal || dist == log_normal){
                printf("using advDir up normal or log_normal dist priors\n");
                if (degree == 1) {
                  prior = prRFreqPoly1Up;
                } else if (degree == 2){
                  prior = prRFreqPoly2Up;
                } else if (degree == 3){
                  prior = prRFreqPoly3Up;
                } else if (degree == 4){
                  prior = prRFreqPoly4Up;
                } else if (degree == 5){
                  prior = prRFreqPoly5Up;
                } else{
                  printf("poly restricted normal/lognormal degree error\n");
                  return;
                }
              } else {
                printf("using advDir up normal_ncv dist priors\n");
                if (degree == 1) {
                  prior = prRFreqPoly1NCVUp;
                } else if (degree == 2){
                  prior = prRFreqPoly2NCVUp;
                } else if (degree == 3){
                  prior = prRFreqPoly3NCVUp;
                } else if (degree == 4){
                  prior = prRFreqPoly4NCVUp;
                } else if (degree == 5){
                  prior = prRFreqPoly5NCVUp;
                } else{
                  printf("poly restricted normal NCV degree error\n");
                  return;
                }
              }

            } else {
              if(dist == normal || dist == log_normal){
                printf("using advDir down normal or log_normal dist priors\n");
                if (degree == 1) {
                  prior = prRFreqPoly1Down;
                } else if (degree == 2){
                  prior = prRFreqPoly2Down;
                } else if (degree == 3){
                  printf("using prRFreqPoly3Down\n");
                  prior = prRFreqPoly3Down;
                } else if (degree == 4){
                  prior = prRFreqPoly4Down;
                } else if (degree == 5){
                  prior = prRFreqPoly5Down;
                } else{
                  printf("poly restricted normal/lognormal degree error\n");
                  return;
                }
              } else {
                printf("using advDir down normal_ncv dist priors\n");
                if (degree == 1) {
                  prior = prRFreqPoly1NCVDown;
                } else if (degree == 2){
                  prior = prRFreqPoly2NCVDown;
                } else if (degree == 3){
                  prior = prRFreqPoly3NCVDown;
                } else if (degree == 4){
                  prior = prRFreqPoly4NCVDown;
                } else if (degree == 5){
                  prior = prRFreqPoly5NCVDown;
                } else{
                  printf("poly restricted normal NCV degree error\n");
                  return;
                }
              }
            }
          }
          break;
        default :
          printf("error with restricted models\n");
          return;

      }
    } else {
      //unrestricted
      switch(model) {

        case hill:
          anal.model = hill;
          if (dist == normal) {
            //normal
            prior = prUFreqHillNormal;
          } else if (dist == normal_ncv){
            //normal NCV
            prior = prUFreqHillNormalNCV;
          } else {
            //lognormal
            prior = prUFreqHillLognormal;
          }
          break;
        case exp_3:
          printf("cannot run unrestricted exponential models\n");
          return;
          //break;
        case exp_5:
          printf("cannot run unrestricted exponential models\n");
          return;
          //break;
        case power:
          anal.model = power;
          if (dist == normal || dist == log_normal){
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
          //if (detectAdvDir){
            if(dist == normal || dist == log_normal){
              printf("prior with normal or lognormal dist\n");
              if (degree == 1) {

                prior = prUFreqPoly1;
              } else if (degree == 2){
                prior = prUFreqPoly2;
              } else if (degree == 3){
                prior = prUFreqPoly3;
              } else if (degree == 4){
                prior = prUFreqPoly4;
              } else if (degree == 5){
                prior = prUFreqPoly5;
              } else{
                printf("poly unrestricted normal/lognormal degree error\n");
                return;
              }
            } else {
              if (degree == 1) {
                prior = prUFreqPoly1NCV;
              } else if (degree == 2){
                prior = prUFreqPoly2NCV;
              } else if (degree == 3){
                prior = prUFreqPoly3NCV;
              } else if (degree == 4){
                prior = prUFreqPoly4NCV;
              } else if (degree == 5){
                prior = prUFreqPoly5NCV;
              } else{
                printf("poly restricted normal NCV degree error\n");
                return;
              }
            }
          //}
          break;

        default :
          printf("error with unrestricted model\n");
          return;

      }
    }
  } else {
  //bayesian
    switch(model) {
       case hill:
         anal.model = hill;
         if (dist == normal || dist == log_normal){
           //normal
           prior = prBayesianHill;
         } else {
           //normal NCV
           prior = prBayesianHillNCV;
         }
         break;
       case exp_3:
         anal.model = exp_3;
         if (dist == normal || dist == log_normal){
           //normal
           prior = prBayesianExp5;
         } else {
           //normal NCV
           prior = prBayesianExp5NCV;
         }
         break;
       case exp_5:
         anal.model = exp_5;
         if (dist == normal || dist == log_normal){
           //normal
           prior = prBayesianExp5;
         } else {
           //normal NCV
           prior = prBayesianExp5NCV;
         }
         break;
       case power:
         anal.model = power;
         if (dist == normal || dist == log_normal){
           //normal
           prior = prBayesianPower;
         } else {
           //normal NCV
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
         if (dist == normal || dist == log_normal){
           //normal
           printf("using Bayesian normal or log_normal dist priors\n");
           if (degree == 1) {
             prior = prBayesianPoly1;
           } else if (degree == 2){
             prior = prBayesianPoly2;
           } else if (degree == 3){
             prior = prBayesianPoly3;
           } else if (degree == 4){
             prior = prBayesianPoly4;
           } else if (degree == 5){
             prior = prBayesianPoly5;
           } else{
             printf("poly restricted normal/lognormal degree error\n");
             return;
           }
         } else {
           //normal NCV
           printf("using Bayesian normal_ncv dist priors\n");
           if (degree == 1) {
             prior = prBayesianPoly1NCV;
           } else if (degree == 2){
             prior = prBayesianPoly2NCV;
           } else if (degree == 3){
             prior = prBayesianPoly3NCV;
           } else if (degree == 4){
             prior = prBayesianPoly4NCV;
           } else if (degree == 5){
             prior = prBayesianPoly5NCV;
           } else{
             printf("poly restricted normal/lognormal degree error\n");
             return;
           }
         }
         break;
    }
  }
  //parms array declared
//  int numParms = sizeof(pr)/sizeof(pr[0])/prCols;
  //double parms[numParms];
  double* parms = new double[numParms];

  //declare analysis
  anal.Y.assign(Y, Y + numDataRows);
  anal.n = numDataRows;
  if(suffStat){
    anal.n_group.assign(N, N + numDataRows);
    anal.sd.assign(SD, SD + numDataRows);
  }
  anal.doses.assign(D, D + numDataRows);
  anal.disttype = dist;
  if (!detectAdvDir){
    anal.isIncreasing = isIncreasing;
  }

  anal.alpha = alpha;
  anal.BMD_type = BMD_type;   //1=absdev, 2 = stddev, 3 = reldev, 4 = pt, 5 = extra, 6 = hybrid_extra, 7 = hybrid_added   from src/include/cmodeldefs.h
  anal.BMR = BMRF;
  anal.samples = 0;   //num MCMC samples
  anal.tail_prob = 0.01;
  anal.suff_stat = suffStat;
  anal.parms = numParms;
  anal.prior_cols = prCols;
  anal.transform_dose = 0;
  anal.prior.assign(prior, prior + anal.prior_cols*anal.parms);
  anal.restricted = restricted;
  anal.detectAdvDir = detectAdvDir;

  printf("prior b4 adj:\n");
  for (int i=0; i<prCols*numParms; i++){
    printf("%.9f\n",anal.prior[i]);
  }


  struct python_continuous_model_result res;
  res.model = anal.model;
  res.nparms = anal.parms;
  res.dist_numE = 100;

  struct BMDS_results BMDSres;
  //set all parms as unbounded initially
  for (int i=0; i<anal.parms; i++){
     BMDSres.bounded.push_back(false);;
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
  if(anal.suff_stat){
    nGOF = anal.n;
  } else {
    //determine number of unique dose groups
    nGOF = 1;
    for (int i=1; i<anal.n; i++){
      int j=0;
      for (j=0; j < i; j++){
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
  double addConst; // = 22.2;
  aod.LL = LL;
  aod.nParms = nParms;
  aod.AIC = AIC;
  aod.addConst = addConst;

  //double llRatio[4];
  //double DF[4];
  //double pVal[4];
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
  printf("priors sent to model code:\n");
  for (int i=0; i<prCols*numParms; i++){
    printf("%.20f\n",anal.prior[i]);
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

 if (anal.suff_stat){
    for (int i=0; i<anal.n; i++){
       printf("%.20f, %.20f, %.20f, %.20f\n",anal.doses[i],anal.n_group[i], anal.Y[i], anal.sd[i]);

    }
  } else {
    for (int i=0; i<anal.n; i++){
       //printf("%.20f, %.20f, %.20f, %.20f\n",anal.doses[i],anal.n_group[i], anal.Y[i], anal.sd[i]);
       printf("%.20f, %.20f\n",anal.doses[i],anal.Y[i]);
    }
  }
  printf("\n\n");

  printf("calling pythonBMDSCont\n");
  pythonBMDSCont(&anal, &res);

  if(detectAdvDir){
    printf("auto adverse direction: %s\n", anal.isIncreasing ? "increasing" : "decreasing");
  }

  printf("\n\n");
  printf("prior after adj by model code:\n");
  for (int i=0; i<prCols*numParms; i++){
    printf("%.20f\n",anal.prior[i]);
  }

  printf("\n\n----------OUTPUT-----------\n");
  printf("tlink BMDSres.validResult = %s\n", res.bmdsRes.validResult ? "valid" : "invalid");
  if (BMDSres.validResult || showResultsOverride){

  printf("\nBenchmark Dose\n");
  printf("max:  %f\n",res.max);
  printf("BMD:  %f\n",res.bmdsRes.BMD);
  printf("Matt's BMD:  %f\n",res.bmd);
  printf("BMDL: %f\n",res.bmdsRes.BMDL);
  printf("BMDU: %f\n",res.bmdsRes.BMDU);
  printf("AIC:  %f\n",res.bmdsRes.AIC);
  printf("LPP: %f\n", res.bmdsRes.BIC_equiv);
  printf("Test 4 P-value: %f\n", aod.TOI.pVal[3]);
  printf("DOF: %f\n", aod.TOI.DF[3]);
  printf("ChiSq: %f\n", res.bmdsRes.chisq);

  printf("\nModel Parameters\n");
  printf("# of parms: %d\n", anal.parms);
  printf("parm, estimate, bounded, std.err., lower conf, upper conf\n");
  for (int i=0; i<anal.parms; i++){
     printf("%d, %.20f, %s, %f, %f, %f\n", i, res.parms[i], res.bmdsRes.bounded[i]? "true" : "false", res.bmdsRes.stdErr[i], res.bmdsRes.lowerConf[i], res.bmdsRes.upperConf[i]);
//     printf("bounded %d = %s\n", i, BMDSres.bounded[i] ? "true" : "false");
  }

  printf("\nGoodness of Fit\n");
  printf("res.gof.n = %d\n",res.gof.n);
  printf("Dose, Size, EstMed, CalcMed, ObsMean, EstSD, CalcSD, ObsSD, SR\n");
  for(int i=0; i<res.gof.n; i++){
    printf("%f, %f, %f, %f, %f, %f, %f, %f, %f\n",res.gof.dose[i],res.gof.size[i],res.gof.estMean[i],res.gof.calcMean[i],res.gof.obsMean[i],res.gof.estSD[i],res.gof.calcSD[i],res.gof.obsSD[i],res.gof.res[i]);
  }
  printf("\nError Bars\n");
  for(int i=0; i<res.gof.n; i++){
    printf("%f, %f\n", res.gof.ebLower[i], res.gof.ebUpper[i]);
  }

  printf("\nLikelihoods of Interest\n");
  for (int i=0; i<5; i++){
    printf("i:%d, LL:%f, nParms:%d, AIC:%f\n",i,res.aod.LL[i],res.aod.nParms[i],res.aod.AIC[i]);
  }
  printf("additive constant:%f\n",res.aod.addConst);

  printf("\nTests of Interest:\n");
  for (int i=0; i<4; i++){
    printf("i:%d, llRatio:%f, DF:%f, pVal:%f\n",i,res.aod.TOI.llRatio[i],res.aod.TOI.DF[i],res.aod.TOI.pVal[i]);
  }

  printf("\nBMD Dist:\n");
  for (int i=0; i<res.dist_numE; i++){
    printf("i:%d, perc:%f, dist:%f\n", i, res.bmd_dist[i+res.dist_numE], res.bmd_dist[i]);
  }

  } else {
    printf("Model was not run\n");
  }


  //debugContAnal(&anal);

}


void runPythonMultitumorAnalysis(){

  enum dich_model model = d_multistage;  //d_hill =1, d_gamma=2,d_logistic=3, d_loglogistic=4,
                                   //d_logprobit=5, d_multistage=6,d_probit=7,
                                   //d_qlinear=8,d_weibull=9
  int modelType = 1;       //1 = frequentist, 2 = bayesian
  bool restricted = true;  //only used for frequentist models
  int BMD_type = 1;        // 1 = extra ; added otherwise
  double BMR = 0.1;
  double alpha = 0.05;

  //data
  std::vector<double> doses1 = {0,50,100,150,200};
  std::vector<double> Y1 = {0,5,30,65,90};
  std::vector<double> n_group1 = {100,100,100,100,100};
  std::vector<double> doses2 = {0,50,100,150,200};
  std::vector<double> Y2 = {5,10,33,67,93};
  std::vector<double> n_group2 = {100,100,100,100,100};
  std::vector<double> doses3 = {0,50,100,150,200};
  std::vector<double> Y3 = {1,68,78,88,98};
  std::vector<double> n_group3 = {100,100,100,100,100};

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

  std::vector<int> n = {5,5,5};
  std::vector<int> degree = {0,0,0};
/////////////////////////////////////////////////
////END USER INPUT
////////////////////////////////////////////////////

  //priors defined columnwise
  int prCols = 5;

  int numDatasets = Y.size();
  //int numDatasets = 1;

  struct python_multitumor_analysis anal;
  anal.ndatasets = numDatasets;
  anal.n = n;
  anal.degree = degree;
  anal.BMR = BMR;
  anal.BMD_type = BMD_type;
  anal.alpha = alpha;
  anal.prior_cols = prCols;

  struct python_multitumor_result res;
  res.ndatasets = numDatasets;

  //create individual models analyses
  std::vector<std::vector<python_dichotomous_analysis>> models;
  int count = 0;



  for (int dataset=0; dataset<numDatasets; dataset++){

    int numDataRows = Y[dataset].size();

    struct python_dichotomous_analysis modAnal;
    modAnal.model = d_multistage;
    printf("model = %d\n",modAnal.model);
    modAnal.BMD_type = BMD_type;
    modAnal.BMR = BMR;
    modAnal.alpha = alpha;
    modAnal.prior_cols = prCols;

    struct python_dichotomous_model_result modRes;
    modRes.model = modAnal.model;
    modRes.dist_numE = 200;


    std::vector<python_dichotomous_analysis> modGroup;
    std::vector<python_dichotomous_model_result> modResGroup;

    //Individual model construction
    //declare analysis
    modAnal.parms = 1 + degree[dataset];
    modAnal.Y = Y[dataset];
    modAnal.n_group = n_group[dataset];
    modAnal.doses = doses[dataset];
    modAnal.n = numDataRows;

    //needs to be changed based on model degree
    if (degree[dataset] == 0){
      //handle autoselect degree
      count = 0;
      //run models from degree = 1 to k-1, where k = number of dose groups
      for (int deg=1; deg<anal.n[dataset]; deg++){
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

  //run MSCombo
  pythonBMDSMultitumor(&anal, &res);


  //individual model results
  for (int dataset=0; dataset<numDatasets; dataset++){
    std::cout<<"dataset:"<<dataset<<std::endl;
    for (int mod=0; mod<anal.nmodels[dataset]; mod++){
        std::cout<<" model:"<<mod<<std::endl;
        printDichoModResult(&anal.models[dataset][mod], &res.models[dataset][mod],true);
    }
  }

  std::cout<<"Selected model Indexes:  ";
  for (auto elem : res.selectedModelIndex) {
        std::cout << elem << ", ";
  }

  std::cout<<"BMD:  "<<res.BMD<<std::endl;
  std::cout<<"BMDL: "<<res.BMDL<<std::endl;
  std::cout<<"BMDU: "<<res.BMDU<<std::endl;
  std::cout<<"slope factor: "<<res.slopeFactor<<std::endl;
  std::cout<<"combined LL: "<<res.combined_LL<<std::endl;
  std::cout<<"combined LL constant: "<<res.combined_LL_const<<std::endl;

}




//Helper routines

void printDichoModResult(struct python_dichotomous_analysis *pyAnal, struct python_dichotomous_model_result *pyRes, bool showResultsOverride){


    printf("tlink bmdsRes.validResult = %s\n", pyRes->bmdsRes.validResult ? "valid" : "invalid");
    if (pyRes->bmdsRes.validResult || showResultsOverride){
       std::cout<<"Valid Result"<<std::endl;
    printf("\nBenchmark Dose\n");
    printf("max: %f\n",pyRes->max);
    printf("BMD: %f\n",pyRes->bmdsRes.BMD);
    printf("BMDL: %f\n",pyRes->bmdsRes.BMDL);
    printf("BMDU: %f\n",pyRes->bmdsRes.BMDU);
    printf("AIC: %f\n",pyRes->bmdsRes.AIC);
    printf("LPP: %f\n", pyRes->bmdsRes.BIC_equiv);
    printf("P-value: %f\n", pyRes->gof.p_value);
    printf("DOF: %f\n", pyRes->gof.df);
    printf("Chi^2: %f\n", pyRes->bmdsRes.chisq);

    printf("\nModel Parameters\n");
    printf("# of parms: %d\n", pyAnal->parms);
    printf("parm, estimate, bounded, std.err., lower conf, upper conf\n");
    for (int i=0; i<pyAnal->parms; i++){
       printf("%d, %.10f, %s, %f, %f, %f\n", i, pyRes->parms[i], pyRes->bmdsRes.bounded[i] ? "true" : "false", pyRes->bmdsRes.stdErr[i], pyRes->bmdsRes.lowerConf[i], pyRes->bmdsRes.upperConf[i] );
    }

    printf("\ncov matrix\n");
    for (int i=0; i<pyAnal->parms*pyAnal->parms; i++){
      printf("%d, %f\n", i, pyRes->cov[i]);
    }

    printf("\nGoodness of Fit\n");
    printf("Dose, EstProb, Expected, Observed, Size, ScaledRes\n");
    for (int i=0; i<pyRes->gof.n; i++){
      printf("%f, %f, %f, %f, %f, %f\n", pyAnal->doses[i], pyRes->gof.expected[i]/pyAnal->n_group[i], pyRes->gof.expected[i], pyAnal->Y[i], pyAnal->n_group[i], pyRes->gof.residual[i]);
    }
    printf("\nError bars\n");
    for (int i=0; i<pyRes->gof.n; i++){
      printf("%f, %f\n", pyRes->gof.ebLower[i], pyRes->gof.ebUpper[i]);
    }

    printf("\nAnalysis of Deviance\n");
    printf("  Model,   LL,    #parms,   deviance,   test DF,  pval\n");
    printf("Full Model,  %f,  %d,  -,  -,  NA\n", pyRes->aod.fullLL, pyRes->aod.nFull);
    printf("Fitted Model,  %f,  %d,  %f,  %d,  %f\n", pyRes->aod.fittedLL, pyRes->aod.nFit, pyRes->aod.devFit, pyRes->aod.dfFit, pyRes->aod.pvFit);
    printf("Reduced Model,  %f,  %d,  %f,  %d,  %f\n", pyRes->aod.redLL, pyRes->aod.nRed, pyRes->aod.devRed, pyRes->aod.dfRed, pyRes->aod.pvRed);

    printf("\nBMD Dist:\n");
    for (int i=0; i<pyRes->dist_numE; i++){
      printf("i:%d, perc:%f, dist:%f\n", i, pyRes->bmd_dist[i+pyRes->dist_numE], pyRes->bmd_dist[i]);
    }
  } else {
     printf("\nModel was not run\n");
  }
}
std::vector<double> getMultitumorPrior(int degree, int prior_cols){

  int numParms = 2;  //2 parameters for multistage cancer G & B
  //initial values for multistage 2
  std::vector<double> prG(prRFreqMultistageCancerG, prRFreqMultistageCancerG + prior_cols);
  std::vector<double> prB(prRFreqMultistageCancerB, prRFreqMultistageCancerB + prior_cols);

  std::vector<double> pr;
  for (int i=0; i<prior_cols; i++){
     pr.push_back(prG[i]);
     for (int j=0; j<degree; j++){
       pr.push_back(prB[i]);
     }
  }
  return pr;

}

