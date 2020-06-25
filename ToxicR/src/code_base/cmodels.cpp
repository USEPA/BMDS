#include "stdafx.h"

#include <fstream>
#include <chrono>
#include <iostream>
#include <stdexcept>

#undef eigen_assert
#define eigen_assert(x) \
  if (!(x)) { throw (std::runtime_error("Eigen library assertion failed.")); }


#ifdef R_COMPILATION
    //necessary things to run in R    
    #include <RcppEigen.h>
    #include <RcppGSL.h>
   // #define _stdcall
#else 
	#include <Eigen/Dense>
	#include <crtdbg.h>
#endif

#include <iomanip>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_randist.h>

#include "bmd_calculate.h"

#include "statmod.h"
#include "log_likelihoods.h"
#include "normal_likelihoods.h"
#include "normalModels.h"
#include "lognormalModels.h"
#include "binomModels.h"
#include "IDPrior.h"

#include "normal_HILL_NC.h"
#include "normal_POWER_NC.h"
#include "normal_POLYNOMIAL_NC.h"
#include "normal_EXP_NC.h"

#include "lognormal_HILL_NC.h"
#include "lognormal_POWER_NC.h"
#include "lognormal_POLYNOMIAL_NC.h"
#include "lognormal_EXP_NC.h"

#include "cBMDstatmod.h"
#include "normalTests.h"
#include "lognormalTests.h"
#include "bmds_entry.h"
#include "bmd_calculate.h"
#include "cModel.h"

#define LOGFILENAME __FUNCTION__ ## ".log"

//#if defined(_DEBUG) || defined(DEBUGLOG)
//#  define DEBUG_LOG(f, s) (f) << (s) << endl;
//#else
//#  define DEBUG_LOG(f, s)
//#endif // _DEBUG


#if defined WIN32 || defined _WINDOWS
#   ifndef R_COMPILATION
#define ALLOC_MODEL_ID(x) SysAllocString(L ## x)
#else
#       define ALLOC_MODEL_ID(x) x
#   endif
#else
//#define ALLOC_MODEL_ID(x) new char*(x)
#define ALLOC_MODEL_ID(x) x
#endif // WIN32

// Base number of model parameters, excluding variance model
#define NPARM_EXP2 2
#define NPARM_EXP3 3
#define NPARM_EXP4 3
#define NPARM_EXP5 4
// Get # parameters from EXP submodel number
#define NPARM_EXP(x)  ((x) < 4 ? (x) : (x) - 1)
#define NPARM_HILL 4
#define NPARM_POLY_BASE 1 // Base # of parms exluding betas according to degree
#define NPARM_POLY(x) (1 + (x))
#define NPARM_POW 3
// Number of variance parameters
#define NPARM_VAR_CV 1  // Constant variance
#define NPARM_VAR_NCV 2 // Modeled variance

#define iBACKGROUND 0 // Array index for background parameter for all models

// Determine EXP model type according to submodel and adverse direction
#define EXP_TYPE(X, Y) ((Y) ? (X) : (X) * 10 + 1)

enum EXPSubmodel_t {eM2 = 2, eM3 = 3, eM4 = 4, eM5 = 5};
static const char* MODELNAME[] = { NULL, NULL, "Exp2", "Exp3", "Exp4", "Exp5", "Hill", "Polynomial", "Power" };

// Starting value and step size for the analysis profile
constexpr auto ANALYSIS_DEFAULT_STEP_SIZE = 0.01; // Step size for generating profile
constexpr auto ANALYSIS_ALPHA_VALUE = 0.01; // Starting value for generating profile

// Starting value and step size to generate CDF values for BMD
constexpr auto BMD_CDF_START = 0.01;
constexpr auto BMD_CDF_STEP = 0.01;

using namespace std;
using namespace Eigen;

// This array stashes the MLEs from EXP submodels 2, 3 & 4 to initialize the
// parm priors for the next submodel.
// ASSUMPTION: the EXP submodels are always run in order!
static double expStashPriors[NPARM_EXP5 + NPARM_VAR_NCV];

// Forward declarations of local functions
int convertDataToSuffStats(BMDSInputData_t *dataIn, int n, bool bLognormal,
  Eigen::MatrixXd& X, Eigen::MatrixXd& Y, double xScale);
void convertLognormalStats(vector <double>& mean, vector <double>& sd, int n);
Eigen::MatrixXd getParmPriors(CModelID_t model, bool bLognormal, bool bConstVar,
  int nparms, const Eigen::MatrixXd& max_parms, int degree, bool isIncreasing,
  int restriction, double yScale);
int setupData(BMDSInputData_t *dataIn, int n, bool *psuff_stat, bool bLognormal,
  Eigen::MatrixXd& X, Eigen::MatrixXd& Y, double *pmaxDose, double *pmaxMean);
Eigen::MatrixXd lognormal_initializer(Eigen::MatrixXd X, Eigen::MatrixXd Y,
  bool suff_stat);
Eigen::MatrixXd normal_initializer(Eigen::MatrixXd X, Eigen::MatrixXd Y,
  bool suff_stat, bool const_var);
//int modelDF(bmd_analysis *fitout, Eigen::MatrixXd priors, vector<bool>& bFixed, bool *zBounded);
void analysisOfDeviance(ContinuousDeviance_t *zOut,
                        Eigen::MatrixXd& Y, Eigen::MatrixXd& X, double llFit,
                        int nEstParms, bool bLognormal, bool bConstVar,
                        bool bSuffStat, int dgCount);

// GLOBAL VARIABLES
// Variables for memory debugging
#ifndef R_COMPILATION
_CrtMemState memStateDiff, memSnapshot1, memSnapshot2, memSnapshot3, memSnapshot4, memSnapshot5;
#endif
/********************** Start of function definitions **********************/

// This is the entry point for all continuous models. ASSUMPTIONS:
// - Caller allocated sufficient space for MAP/MLE values
int _stdcall run_cmodel2(CModelID_t *p_model, BMD_C_ANAL *returnV, BMDSInputType_t *p_inputType,
  BMDSInputData_t *dataIn, PRIOR *priorsIn, BMDS_C_Options_t *options, int *p_n,
  bmd_analysis *pAnalysis) {
  int exitStatus = -1; // Default to error status
  gsl_set_error_handler_off();
  CModelID_t model = *p_model;
  BMDSInputType_t inputType = *p_inputType;
  int n = *p_n; // # data rows (can change if non-ss data converted to ss)
  int nOrig = *p_n; // # data rows (original count)
  cModel *pModel = NULL;
  int degree = options->degree; // Specifies poly degree OR EXP2, EXP3, EXP4 or EXP5
  bool bLognormal = options->bLognormal;
  bool bConstVar;
  bool bExpModel = false;
  contbmd bmrType = (contbmd)options->bmrType;
  double cl_alpha = options->alpha; // Confidence limit value from user
  double a_alpha = ANALYSIS_ALPHA_VALUE; // Starting value for profile
  double bmrf = options->bmr;
  double tail_prob = options->tailProb;
  double step = ANALYSIS_DEFAULT_STEP_SIZE; // LCO: Add step to BMDS_C_Options_t
  double background = options->background;
  int nparms;
  int adverseDir = options->adverseDirection;
  int restriction = options->restriction;
  bool isIncreasing;
  bool bUserParmInit = options->bUserParmInit;
  string myName(MODELNAME[model]);
  //const char *pcModelName = MODELNAME[model];

//#define DEBUG_HEAP
#ifdef DEBUG_HEAP
  // Enable heap memory debugging - only useful when compiled with -D_DEBUG
  int crtDebugFlags;
  crtDebugFlags = _CrtSetDbgFlag(_CRTDBG_REPORT_FLAG);
  // By default, memory is checked every 1024 allocation calls
  // - Use _CRTDBG_CHECK_ALWAYS_DF to check every call (REALLY SLOW!!!!)
  // - Alternatively, clear the upper 16 bits and "OR" in the desired frequency:
  //   e.g., ( crtDebugFlags & 0x0000FFFF) | _CRTDBG_CHECK_EVERY_16_DF
  // crtDebugFlags = crtDebugFlags | _CRTDBG_CHECK_ALWAYS_DF
  crtDebugFlags = crtDebugFlags | _CRTDBG_CHECK_CRT_DF
    //| _CRTDBG_LEAK_CHECK_DF
    | _CRTDBG_DELAY_FREE_MEM_DF;
  crtDebugFlags = (crtDebugFlags & 0x0000FFFF)
    //| 0x04000000; // every 1024 calls
    //| 0x10000000; // every 4096 calls
    //| 0x40000000; // every 16384 calls
    | 0x80000000; // every 32768 calls

    _CrtSetDbgFlag(crtDebugFlags);

  // Send all reports to STDOUT
  _CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_FILE);
  _CrtSetReportFile(_CRT_WARN, _CRTDBG_FILE_STDOUT);
  _CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_FILE);
  _CrtSetReportFile(_CRT_ERROR, _CRTDBG_FILE_STDOUT);
  _CrtSetReportMode(_CRT_ASSERT, _CRTDBG_MODE_FILE);
  _CrtSetReportFile(_CRT_ASSERT, _CRTDBG_FILE_STDOUT);

#endif // DEBUG_HEAP

  // Initialize return values
  returnV->model_id = ALLOC_MODEL_ID("1.1");
  returnV->bAdverseUp = 0;
  returnV->MAP = 0;
  returnV->BMD = -9999;
  returnV->BMDL = -9999;
  returnV->BMDU = -9999;
  returnV->nparms = 0;

  if (bLognormal) bConstVar = true;
  else bConstVar = (options->varType == eConstant);

  switch (model) {
  case eExp2:
  case eExp3:
  case eExp4:
  case eExp5:
    bExpModel = true;
    nparms = NPARM_EXP5;
    break;
  case eHill:
    nparms = NPARM_HILL;
    break;
  case ePoly:
    nparms = NPARM_POLY_BASE + degree;
    myName.append('-' + to_string(degree));
    break;
  case ePow:
    nparms = NPARM_POW;
    break;
  default:
    break;
  } // end switch
  //int nparms_actual = NPARM_EXP(submodel);

  if (bConstVar) nparms += NPARM_VAR_CV;
  else nparms += NPARM_VAR_NCV;

  bool suff_stat = (inputType == eCont_4); // suff_stat means summarized data

#if defined(_DEBUG) || defined(DEBUGLOG)
  ofstream file;
  string fname(__FUNCTION__);
  fname.append('.' + myName + ".log");
  file.open(fname/*, fstream::app*/);
  time_t now = time(0);
  struct tm localtm;
  char str[64]; memset(str, '\0', sizeof(str));
#ifndef R_COMPILATION
  localtime_s(&localtm, &now);
  asctime_s(str, sizeof(str), &localtm);
#endif
  file << "Running model " << myName << " at " << str;
  file << "\tbLognormal= " << bLognormal << " bConstVar= " << bConstVar
    << " suff_stat= " << suff_stat << " degree= " << degree << "\n";
  file << "\tbmrType= " << bmrType << " bmrf= " << bmrf << " alpha= "
    << cl_alpha << ", background= " << background << " restriction= " << restriction
    << "\n\tstep size= " << ANALYSIS_DEFAULT_STEP_SIZE
    << " tail_prob= " << tail_prob << endl;
  file << "bUserParmInit= " << bUserParmInit << ", nparms= " << nparms << endl;
#endif // _DEBUG

  // Always allocate dose-response data matrices. Lengths can change if
  // individual responses are converted to sufficient stats
  Eigen::MatrixXd  X(n, 1); Eigen::MatrixXd  Y(n, 3);
  // Scale doses - always
  // Scale means - not currently done
  double maxDose = 0, yScale = 1;
  n = setupData(dataIn, n, &suff_stat, bLognormal, X, Y, &maxDose, &yScale);
#if defined(_DEBUG) || defined(DEBUGLOG)
  //file << std::scientific << setprecision(17);
  file << "Input data with doses scaled by: " << maxDose << endl;
  if (suff_stat) {
    file << " Dose\tResp\tStdDev\tSize" << endl;
    for (int i = 0; i < n; i++) {
      file << " " << X(i, 0)
        << '\t' << Y(i, 0)
        << '\t' << Y(i, 1)
        << '\t' << Y(i, 2)
        << '\n';
    } // for
  }
  else {
    file << " Dose\tResp" << endl;
    for (int i = 0; i < n; i++) {
      file << " " << X(i, 0)
        << '\t' << Y(i, 0)
        << '\n';
    } // for
  }
  file << endl;
#endif // _DEBUG

  Eigen::MatrixXd max_parms;
  //double bParm, *pbParm; // Only used for EXP model
  //pbParm = bExpModel ? &bParm : NULL;
  if (bLognormal) max_parms = lognormal_initializer(X, Y, suff_stat);
  else max_parms = normal_initializer(X, Y, suff_stat, bConstVar);
  //if (pbParm != NULL) *pbParm *= maxDose; // Scale by dose

  // direction of adversity according to the regression
  returnV->bAdverseUp = isIncreasing = (max_parms(1, 0) > 0);
  // Override direction of adversity if specified by the user
  if (adverseDir == 1) isIncreasing = true;
  else if (adverseDir == -1) isIncreasing = false;
  // Handle automatic up/down restriction for the poly model
  if (restriction == BMDS_BLANK_VALUE) restriction = isIncreasing ? 1 : -1;

  // Only the exponential model is valid for lognormal and adverse-down data
  if ( bLognormal && !isIncreasing && !bExpModel) {
    DEBUG_LOG(file, "Error - Only the exponential models are valid for Lognormal and adverse down data");
    return exitStatus;
  }

#if defined(_DEBUG) || defined(DEBUGLOG)
  file << "Parameter initializer vector:" << endl;
  file << "initializer #rows= " << max_parms.rows() << ", #cols= " << max_parms.cols() << endl;
  for (int i = 0; i < max_parms.rows(); i++) {
    file << max_parms(i, 0) << endl;
  }
  file << "\n" << endl;
#endif // _DEBUG

  // If bUserParmInit is true, use the passed-in prior/initial values;
  // otherwise, use automatically calculated values.
  Eigen::MatrixXd initPriors(nparms, NUM_PRIOR_COLS);
  if (bUserParmInit) {
    for (int i = 0; i < nparms; i++) {
      initPriors(i, 0) = priorsIn[i].type;
      initPriors(i, 1) = priorsIn[i].initalValue;
      initPriors(i, 2) = priorsIn[i].stdDev;
      initPriors(i, 3) = priorsIn[i].minValue;
      initPriors(i, 4) = priorsIn[i].maxValue;
    } // for
  }
  else {
    initPriors = getParmPriors(model, bLognormal, bConstVar,
      nparms, max_parms, degree,
      isIncreasing, restriction, yScale);
  }
#if defined(_DEBUG) || defined(DEBUGLOG)
  file << "Parameter initial values / priors:" << endl;
  file << "Type\tInitVal\tStDev\tMinVal\tMaxVal" << endl;
  for (int i = 0; i < nparms; i++) {
    file << initPriors(i, 0) << '\t' << initPriors(i, 1) << '\t' << initPriors(i, 2)
         << '\t' << initPriors(i, 3) << '\t' << initPriors(i, 4) << endl;
  }
  file << endl;
#endif // _DEBUG

  // Fixed parameters are specified w/r/t priors and not the actual parm list
  std::vector<bool> fixedB(nparms);
  std::vector<double> fixedV(nparms);
  for (int i = 0; i < nparms; i++) {
    fixedB[i] = false;
    fixedV[i] = 0.0;
  }
  if (background != BMDS_BLANK_VALUE) {
    //file << "Fixing parameter# " << iBACKGROUND << " to " << background << endl;
    initPriors(0, 1) = fixedV[iBACKGROUND] = background;
    fixedV[iBACKGROUND] = background;
    fixedB[iBACKGROUND] = true;
  }

  if (!isIncreasing && bmrType == (int)BMRType_t::eRelativeDev) bmrf = 1 - bmrf;
  IDcontinuousPrior zPriors(initPriors);
  int TYPE = 0; // Needed by exponential bmd classes
  bmd_analysis a;

#ifdef DEBUG_HEAP
  cout << "***** Calling checkMemory @ line " << __LINE__ << endl;
  _CrtCheckMemory();
  _CrtMemCheckpoint(&memSnapshot1); // Memory debugging
#endif // DEBUG_HEAP

  try {
    switch (model) {
    case eExp2:
    case eExp3:
    case eExp4:
    case eExp5:
      TYPE = EXP_TYPE(model, isIncreasing);
      if (bLognormal) {
        lognormalEXPONENTIAL_BMD_NC likelihood(Y, X, suff_stat, TYPE);
        a = bmd_analysis_CNC<lognormalEXPONENTIAL_BMD_NC, IDcontinuousPrior>
          (likelihood, zPriors, fixedB, fixedV, bmrType, bmrf, tail_prob, isIncreasing, a_alpha, step);
      }
      else {
        normalEXPONENTIAL_BMD_NC likelihood(Y, X, suff_stat, bConstVar, TYPE);
        a = bmd_analysis_CNC<normalEXPONENTIAL_BMD_NC, IDcontinuousPrior>
          (likelihood, zPriors, fixedB, fixedV, bmrType, bmrf, tail_prob, isIncreasing, a_alpha, step);
      } // endif (bLognormal)
      break;
    case eHill:
      if (bLognormal) {
        lognormalHILL_BMD_NC likelihood(Y, X, suff_stat, 0);
        a = bmd_analysis_CNC<lognormalHILL_BMD_NC, IDcontinuousPrior>
          (likelihood, zPriors, fixedB, fixedV, bmrType, bmrf, tail_prob, isIncreasing, a_alpha, step);
      }
      else {
        normalHILL_BMD_NC likelihood(Y, X, suff_stat, bConstVar, 0);
        a = bmd_analysis_CNC<normalHILL_BMD_NC, IDcontinuousPrior>
          (likelihood, zPriors, fixedB, fixedV, bmrType, bmrf, tail_prob, isIncreasing, a_alpha, step);
      } // endif (bLognormal)
      break;
    case ePoly:
      if (bLognormal) {
        lognormalPOLYNOMIAL_BMD_NC likelihood(Y, X, suff_stat, degree);
        a = bmd_analysis_CNC<lognormalPOLYNOMIAL_BMD_NC, IDcontinuousPrior>
          (likelihood, zPriors, fixedB, fixedV, bmrType, bmrf, tail_prob, isIncreasing, a_alpha, step);
      }
      else {
        normalPOLYNOMIAL_BMD_NC likelihood(Y, X, suff_stat, bConstVar, degree);
        a = bmd_analysis_CNC<normalPOLYNOMIAL_BMD_NC, IDcontinuousPrior>
          (likelihood, zPriors, fixedB, fixedV, bmrType, bmrf, tail_prob, isIncreasing, a_alpha, step);
      } // endif (bLognormal)
      break;
    case ePow:
      if (bLognormal) {
        lognormalPOWER_BMD_NC likelihood(Y, X, suff_stat, 0);
        a = bmd_analysis_CNC<lognormalPOWER_BMD_NC, IDcontinuousPrior>
          (likelihood, zPriors, fixedB, fixedV, bmrType, bmrf, tail_prob, isIncreasing, a_alpha, step);
      }
      else {
        normalPOWER_BMD_NC likelihood(Y, X, suff_stat, bConstVar, 0);
        a = bmd_analysis_CNC<normalPOWER_BMD_NC, IDcontinuousPrior>
          (likelihood, zPriors, fixedB, fixedV, bmrType, bmrf, tail_prob, isIncreasing, a_alpha, step);
      } // endif (bLognormal)
      break;
    default:
      break;
    } // end switch
  } // try
  catch (std::exception &e) {
    DEBUG_LOG(file, "Error - analysis failed with an unexpected exception:");
    DEBUG_LOG(file, e.what());
    return exitStatus;
  } // catch
  catch (...) {
    DEBUG_LOG(file, "Error - analysis failed with an unknown exception:");
    return exitStatus;
  } // catch

  // Calculate the additive constant value that is omitted by BMDS 2.7 but
  // included in these likelihoods. Return the constant to simplify
  // comparisons with BMDS 2.7.
  double additiveConstant;
  additiveConstant = (suff_stat) ? -Y.col(2).sum() * log(2 * M_PI) / 2
                                 : -X.rows() * log(2 * M_PI) / 2;
  if (bLognormal) {
    // The additive constant has an additional component for lognormal
    if (suff_stat) additiveConstant -= (Y.col(2).array()*Y.col(0).array()).sum();
    else           additiveConstant -= Y.col(0).sum();
  } // if (bLognormal)
#if defined(_DEBUG) || defined(DEBUGLOG)
  file << "Model Results BEFORE re-scaling:" << endl;
  file << " MAP (-LL) = " << a.MAP << '\n';
  file << " BMD= " << a.MAP_BMD << " BMDL= " << a.BMD_CDF.inv(cl_alpha) << " BMDU= "
    << a.BMD_CDF.inv(1.0 - cl_alpha) << '\n';
  file << "MLE Parameter values\n";
  for (int i = 0; i < nparms; i++) {
    file << "Parm" << i << '\t' << a.MAP_ESTIMATE(i, 0) << endl;
  }
  file << "\n" << endl;
  file << " Covariate matrix\n";
  file << a.COV << endl;
#endif // _DEBUG

  // Populate return values to the caller. Some values are re-scaled now;
  // model-specific adjustments occur further down.
  returnV->MAP = -a.MAP; // a.MAP is the negative log-likelihood
  returnV->ll_const = additiveConstant;
  // Multiply by maaxDose to rescale BMDs
  returnV->BMD = a.MAP_BMD * maxDose;
  returnV->BMDL = a.BMD_CDF.inv(cl_alpha) * maxDose;
  returnV->BMDU = a.BMD_CDF.inv(1.0 - cl_alpha) * maxDose;
  returnV->nparms = nparms;
  for (int i = 0; i < nparms; i++) {
    returnV->PARMS[i] = a.MAP_ESTIMATE(i, 0);
  }
  // Re-scale the variance parameter.
  //returnV->PARMS[nparms - 1] += log(yScale*yScale);
  //if (!bConstVar) {
  //  returnV->PARMS[nparms - 1] -= returnV->PARMS[nparms - 2] * log(yScale);
  //}

  // Populate CDF values to return
  double myVal = BMD_CDF_START;
  for (int i = 0; i < returnV->nCDF; i++) {
    returnV->aCDF[i] = a.BMD_CDF.inv(myVal) * maxDose;
    myVal += BMD_CDF_STEP;
  }

  // Deal with potential exceptions related to goodness-of-fit calculations
  try {
    // We do two things here:
    // 1. For EXP models, stash certain MLEs to init priors for subsequent EXP submodels
    // 2. For all models, re-scale parameters as necessary

    switch (model) {
    case eExp2:
      expStashPriors[0] = a.MAP_ESTIMATE(0, 0); // a
      expStashPriors[1] = a.MAP_ESTIMATE(1, 0); // b
      expStashPriors[nparms-1] = a.MAP_ESTIMATE(nparms-1, 0); // log-alpha
      // Re-scale parameters being returned
      returnV->PARMS[0] *= yScale;  // a
      returnV->PARMS[1] /= maxDose; // b
      // create model object for GoF calculations
      pModel = new bmds_exp2(bConstVar, bLognormal);
      break;
    case eExp3:
      //expStashPriors[1] = a.MAP_ESTIMATE(1, 0); // b
      expStashPriors[3] = a.MAP_ESTIMATE(3, 0); // d
      //expStashPriors[nparms-1] = a.MAP_ESTIMATE(nparms - 1, 0); // ln-alpha
      // Re-scale parameters being returned
      returnV->PARMS[0] *= yScale;  // a
      returnV->PARMS[1] /= maxDose; // b
      // create model object for GoF calculations
      pModel = new bmds_exp3(bConstVar, bLognormal);
      break;
    case eExp4:
      expStashPriors[0] = a.MAP_ESTIMATE(0, 0); // a
      //expStashPriors[1] = a.MAP_ESTIMATE(1, 0); // b
      expStashPriors[2] = a.MAP_ESTIMATE(2, 0); // c
      // Re-scale parameters being returned
      returnV->PARMS[0] *= yScale;  // a
      returnV->PARMS[1] /= maxDose;               // b
      returnV->PARMS[2] = exp(returnV->PARMS[2]); // c
      // create model object for GoF calculations
      pModel = new bmds_exp4(bConstVar, bLognormal);
      break;
    case eExp5:
      // Re-scale parameters being returned
      returnV->PARMS[0] *= yScale;  // a
      returnV->PARMS[1] /= maxDose; // b
      returnV->PARMS[2] = exp(returnV->PARMS[2]); // c
      // create model object for GoF calculations
      pModel = new bmds_exp5(bConstVar, bLognormal);
      break;
    case eHill:
      // Re-scale parameters being returned
      returnV->PARMS[2] *= maxDose; // k (slope)
      // create model object for GoF calculations
      pModel = new bmds_hill(bConstVar, bLognormal);
      break;
    case ePoly:
      // Re-scale parameters being returned - only the betas
      for (int i = 1; i <= degree; i++) {
        returnV->PARMS[i] /= pow(maxDose, i);
      }
      // create model object for GoF calculations
      pModel = new bmds_poly(degree, bConstVar, bLognormal);
      break;
    case ePow:
      // Re-scale parameters being returned
      returnV->PARMS[1] /= pow(maxDose, returnV->PARMS[2]); // slope
      // create model object for GoF calculations
      pModel = new bmds_power(bConstVar, bLognormal);
      break;
    default:
      break;
    } // end switch

#if defined(_DEBUG) || defined(DEBUGLOG)
    file << "Model Results:" << endl;
    file << " BMD= " << returnV->BMD << " BMDL= " << returnV->BMDL <<
      " BMDU= " << returnV->BMDU << '\n';
    file << " bAdverseUp= " << returnV->bAdverseUp << '\n';
    for (int i = 0; i < returnV->nparms; i++) {
      file << "Parm" << i << '\t' << returnV->PARMS[i] << endl;
    }
    file << endl;
#endif // _DEBUG

#ifdef DEBUG_HEAP
    cout << "***** Calling checkMemory @ line " << __LINE__ << endl;
    _CrtCheckMemory();
    _CrtMemCheckpoint(&memSnapshot1); // Memory debugging
#endif // DEBUG_HEAP

    // Perform goodness-of-fit related calculations
    pModel->setParms(returnV->PARMS, isIncreasing);
    // Convert lnalpha to alpha, if appropriate.
    if (bConstVar && !bExpModel) {
      returnV->PARMS[nparms - 1] = exp(returnV->PARMS[nparms - 1]);
    }
#ifdef DEBUG_HEAP
    cout << "***** Calling checkMemory @ line " << __LINE__ << endl;
    _CrtCheckMemory();
    //_CrtMemCheckpoint(&memSnapshot2); // Memory debugging
    //cout << "***** memDiff 1->2..." << endl;
    //if (_CrtMemDifference(&memStateDiff, &memSnapshot1, &memSnapshot2)) {
    //  _CrtMemDumpStatistics(&memStateDiff);
    //}
#endif // DEBUG_HEAP

    pModel->setPriors(initPriors);
#ifdef DEBUG_HEAP
    cout << "***** Calling checkMemory @ line " << __LINE__ << endl;
    _CrtCheckMemory();
    //_CrtMemCheckpoint(&memSnapshot3); // Memory debugging
    //cout << "***** memDiff 2->3..." << endl;
    //if (_CrtMemDifference(&memStateDiff, &memSnapshot2, &memSnapshot3)) {
    //  _CrtMemDumpStatistics(&memStateDiff);
    //}
#endif // DEBUG_HEAP
    int estParmCount; // Number of model parameters that are estimated
    estParmCount = pModel->modelDF(a.MAP_ESTIMATE, fixedB, returnV->boundedParms);
    int dgCount = 0; // Number of dose groups
    dgCount = pModel->calcGoF(returnV->gofRow, inputType, dataIn, nOrig, estParmCount, X, Y);
    // Calculate AIC value for frequentist runs
    double aic;
    // returnV->MAP has the LL not MAP (where MAP is -LL)
    // BMDS 2.7 excludes the additive constant
    aic = -2 * (returnV->MAP - estParmCount);
    returnV->AIC = aic;
    analysisOfDeviance(&(returnV->deviance), Y, X, returnV->MAP, estParmCount,
      bLognormal, bConstVar, suff_stat, dgCount);

#if defined(_DEBUG) || defined(DEBUGLOG)
    file << "Est. parm count = " << estParmCount << " LL additive const = "
      << additiveConstant << endl;
    file << "BMDS 3.x: LL= " << returnV->MAP << ", AIC= " << aic << endl;
    file << "BMDS 2.7: LL= " << returnV->MAP - additiveConstant << ", AIC= "
      << -2 * (returnV->MAP - additiveConstant - estParmCount) << endl;
#endif // _DEBUG...

  } // end try
  catch (std::exception &e) {
    cerr << "Error - analysis of deviance failed with an unexpected exception:\n";
    cerr << e.what();
    flush(cerr);
    DEBUG_LOG(file, "Error - analysis of deviance failed with an unexpected exception:");
    DEBUG_LOG(file, e.what());
    return exitStatus;
  } // catch
  catch (...) {
    cerr << "Error - analysis of deviance failed with an unknown exception:\n";
    flush(cerr);
    DEBUG_LOG(file, "Error - analysis of deviance failed with an unknown exception:");
    return exitStatus;
  } // catch

#if defined(_DEBUG) || defined(DEBUGLOG)
  file << "Tests of interest:" << endl;
  file << "Test#" << "\tdeviance" << "\tdf " << "\tp-value" << endl;
  for (int i = 0; i < NUM_TESTS_OF_INTEREST; i++) {
    file << returnV->deviance.testRows[i].testNumber
      << '\t' << returnV->deviance.testRows[i].deviance
      << "\t\t" << returnV->deviance.testRows[i].df
      << '\t' << returnV->deviance.testRows[i].pvalue
      << endl;
  }

  file << "Exiting run_cmodel()\n" << endl;
  file.close(); // Close the log file
#endif // _DEBUG

  if (pAnalysis != NULL) *pAnalysis = a;

  exitStatus = 0; // Success
  return exitStatus;
}

int _stdcall run_cmodel(CModelID_t *p_model, BMD_C_ANAL *returnV, BMDSInputType_t *p_inputType,
  BMDSInputData_t *dataIn, PRIOR *priorsIn, BMDS_C_Options_t *options, int *p_n) {
  int retval = -1;
  bmd_analysis a;

  retval = run_cmodel2(p_model, returnV, p_inputType, dataIn,
    priorsIn, options, p_n, &a);

  return (retval);
}

// setupData assumes that X and Y matrices were allocated by the caller. It also:
// - Scales dose values, which improves numerical robustness
// - Converts invidual d-r data into sufficient statistics, which also
//   seems to improve convergence.
// - Returns the new data row count
int setupData(BMDSInputData_t *dataIn, int n, bool *psuff_stat, bool bLognormal,
  Eigen::MatrixXd& X, Eigen::MatrixXd& Y, double *pmaxDose, double *pyScale)
{
  int rowCount = n; // Number of input data rows to be returned

  // Find max dose to scale doses: 0 < dose[i] <= 1
  // Also find min dose, used later to get the response scale factor

  double dmax = 0, dMin = DBL_MAX, yScale=1;
  for (int i = 0; i < n; i++) {
    if (dataIn[i].dose > dmax) dmax = dataIn[i].dose;
    else if (dataIn[i].dose < dMin) dMin = dataIn[i].dose;
  }

  if (*psuff_stat) { // Summarized dose-response data
    vector <double> mean(n);
    vector <double> sd(n);
    for (int i = 0; i < n; i++) {
      mean[i] = dataIn[i].response;
      sd[i] = dataIn[i].col4;
    } // for
    if (bLognormal) convertLognormalStats(mean, sd, n);
    for (int i = 0; i < n; i++) {
      X(i, 0) = dataIn[i].dose / dmax;
      Y(i, 0) = mean[i];
      Y(i, 1) = sd[i];
      Y(i, 2) = dataIn[i].groupSize;
    } // for
  }
  else { // individual dose-response data
    rowCount = convertDataToSuffStats(dataIn, n, bLognormal, X, Y, dmax);
    *psuff_stat = true;
  } // endif (suff_stat)

  *pmaxDose = dmax;
  *pyScale = yScale;
  return (rowCount);
}

Eigen::MatrixXd normal_initializer(Eigen::MatrixXd X, Eigen::MatrixXd Y,
  bool suff_stat, bool const_var)
{

  //DEBUG_OPEN_LOG("bmds.log", file);
  std::vector<double> fixed1(3); for (int i = 0; i < fixed1.size(); i++) { fixed1[i] = 0.0; }
  std::vector<double> fixed2(4); for (int i = 0; i < fixed2.size(); i++) { fixed2[i] = 0.0; }

  std::vector<bool> isfixed1(3); for (int i = 0; i < isfixed1.size(); i++) { isfixed1[i] = false; }
  std::vector<bool> isfixed2(4); for (int i = 0; i < isfixed2.size(); i++) { isfixed2[i] = false; }

  /////////////////////////////////////////////////////////////////
  // Initialize using a simple linear regression
  // The simple linear regression estimates can then 
  // be mapped to the given parameters as starting values
  /////////////////////////////////////////////////////////////////
  MatrixXd ipriors(3, 5);
  ipriors << 0, 37, 1, -1e8, 1e8,
    0, 0, 1, -1e8, 1e8,
    0, 0, 1, -1e8, 1e8;

  IDcontinuousPrior 			 initPrior(ipriors);
  normalPOLYNOMIAL_BMD_NC		 startingV(Y, X, suff_stat, true, 1);
  normalPOLYNOMIAL_BMD_NC		 startingV2(Y, X, suff_stat, false, 1);

  cBMDModel<normalPOLYNOMIAL_BMD_NC, IDcontinuousPrior> model(startingV, initPrior,
    isfixed1, fixed1, true);
  optimizationResult initoR = findMAP<normalPOLYNOMIAL_BMD_NC, IDcontinuousPrior>(&model);
  //DEBUG_LOG(file, "initoR: result= " << initoR.result << ", functionV= " << initoR.functionV
  //  << ", max_parms:\n" << initoR.max_parms);

  MatrixXd ipriors2(4, 5);
  ipriors2 << 0, initoR.max_parms(0, 0), 1, 0, 1e8,
    0, initoR.max_parms(1, 0), 1, -100, 100,
    0, 0, 1, -100, 100,
    0, initoR.max_parms(2, 0), 1, -100, 100;

  IDcontinuousPrior 			 initPrior2(ipriors2);
  cBMDModel<normalPOLYNOMIAL_BMD_NC, IDcontinuousPrior> model2(startingV2, initPrior2,
    isfixed2, fixed2, true);

  optimizationResult initoR2 = findMAP<normalPOLYNOMIAL_BMD_NC, IDcontinuousPrior>(&model2);
  //DEBUG_LOG(file, "initoR2: result= " << initoR2.result << ", functionV= " << initoR2.functionV
  //  << ", max_parms:\n" << initoR2.max_parms);

  //// return the correct matrix
  //DEBUG_CLOSE_LOG(file);
  return const_var ? initoR.max_parms : initoR2.max_parms;
}

/********************************************************************************/
/* This routine converts individual dose-response data to sufficient statistics */
/********************************************************************************/
// NOTE: Remember to reset the suff_stat flag in the caller!

int convertDataToSuffStats(BMDSInputData_t *dataIn, int n, bool bLognormal,
                           Eigen::MatrixXd& X, Eigen::MatrixXd& Y, double xScale)
{
  int ret = 0;

  // Find max dose to scale doses: 0 < dose[i] <= 1
  // Also find min dose, used later to get the response scale factor

  // *** Assume that doses are sorted starting with dMin ***
  // First pass: Convert responses to lognormal if applicable
  vector <double> response(n);
  if (bLognormal) {
    for (int i = 0; i < n; i++) response[i] = log(dataIn[i].response);
  }
  else {
    for (int i = 0; i < n; i++) response[i] = dataIn[i].response;
  } // end if (bLognormal)

  // Third pass: summarize the data
  // Setup and prime the loop variables
  vector<double> vDose, vMean, vSD, vN;
  int rowCount = 1; // # rows in summarized data set
  int gSize = 1; // Number of observations in a dose group
  double dose = dataIn[0].dose;
  double sum1 = response[0]; // sum of responses in group
  double sum2 = response[0] * response[0]; // sum squares of responses
  // NOTE: Loop starts at the 2nd element!!!!
  for (int i = 1; i < n; i++) {
    if (dataIn[i].dose != dose) { // New dose group
      // Wrap up the current group
      vDose.push_back(dose);
      vMean.push_back(sum1 / gSize);
      vN.push_back(gSize);
      vSD.push_back(sqrt((sum2 - sum1 * sum1 / gSize) / (gSize - 1)));
      // Start a new dose group
      dose = dataIn[i].dose;
      ++rowCount;
      gSize = 0;
      sum1 = sum2 = 0;
    } // if
    ++gSize;
    sum1 += response[i];
    sum2 += response[i] * response[i];
  } // for
  // Wrap up the last dose group
  vDose.push_back(dose);
  vMean.push_back(sum1 / gSize);
  vN.push_back(gSize);
  vSD.push_back(sqrt((sum2 - sum1 * sum1 / gSize) / (gSize - 1)));

  // Resize the input data vectors to match the summarized values
  X.resize(rowCount, 1);
  Y.resize(rowCount, 3);

  // save the new data values
  for (int i = 0; i < rowCount; i++) {
    X(i, 0) = vDose[i] / xScale;
    Y(i, 0) = vMean[i]; Y(i, 1) = vSD[i]; Y(i, 2) = vN[i];
  }

  //std::cout << "Summarized dose-response data" << endl;
  //std::cout << "x vecor:\n" << x2 << endl;
  //std::cout << "y matrix:\n" << y2 << endl;
  //flush(std::cout);
  ret = rowCount;
  return(ret);
}

/********************************************************************************/
/* This routine converts individual dose-response data to sufficient statistics */
/* ALTERNATE VERSION */
/********************************************************************************/
int convertDataToSuffStatsV0(Eigen::MatrixXd X, Eigen::MatrixXd Y, double *pb_parm)
{
  int ret = 0;

  MatrixXd x2, y2;
    int nrows = X.rows();
    vector<double> doses, means, stdevs, gCounts;
    double dose;
    double stdev;
    double sum1; // Sum of responses while processing a dose group
    double sum2; // Sum of squares of responses while processing a dose grp
    int dCount = 0; // Count of dose groups
    int gSize = 1; // Number of observations in a dose group

    // Calculate stats in case first group has only one dose
    dose = X(0, 0);
    sum1 = Y(0, 0);
    sum2 = Y(0, 0) * Y(0, 0);
    // NOTE: Loop starts at the 2nd element!!!!
    for (int i = 1; i < nrows; i++) {
      if (X(i, 0) != dose) {
        ++dCount;
        doses.push_back(dose);
        means.push_back(sum1 / gSize);
        gCounts.push_back(gSize);
        stdevs.push_back(sqrt((sum2 - sum1 * sum1 / gSize) / (gSize - 1)));
        gSize = 1;
        dose = X(i, 0);
        sum1 = Y(i, 0);
        sum2 = Y(i, 0) * Y(i, 0);
      }
      else {
        ++gSize;
        sum1 += Y(i, 0);
        sum2 += Y(i, 0) * Y(i, 0);
      } // if
    } // for
    // Process the last dose group
    ++dCount;
    doses.push_back(dose);
    means.push_back(sum1 / gSize);
    gCounts.push_back(gSize);
    stdevs.push_back(sqrt((sum2 - sum1 * sum1 / gSize) / (gSize - 1)));
    //Map<VectorXd> myX(&doses[0], dCount);
    MatrixXd myX(dCount, 1), myY(dCount, 3);
    for (int i = 0; i < dCount; i++) {
      myX(i, 0) = doses[i];
      myY(i, 0) = means[i];  myY(i, 1) = stdevs[i]; myY(i, 2) = gCounts[i];
    }
    x2 = myX;
    y2 = myY;
  
  //std::cout << "Summarized dose-response data" << endl;
  //std::cout << "x vecor:\n" << x2 << endl;
  //std::cout << "y matrix:\n" << y2 << endl;
  //flush(std::cout);


  return(ret);
}

Eigen::MatrixXd lognormal_initializer(Eigen::MatrixXd X, Eigen::MatrixXd Y,
  bool suff_stat)
{
  std::vector<double> fixed1(3); for (int i = 0; i < fixed1.size(); i++) { fixed1[i] = 0.0; }
  std::vector<bool> isfixed1(3); for (int i = 0; i < isfixed1.size(); i++) { isfixed1[i] = false; }
  /////////////////////////////////////////////////////////////////
  // Initialize using a simple linear regression
  // The simple linear regression estimates can then 
  // be mapped to the given parameters as starting values
  /////////////////////////////////////////////////////////////////
  MatrixXd ipriors(3, 5);
  ipriors << 0, 1, 1, 0, 1e8,
    0, 0, 1, -1000, 1000,
    0, 0, 1, -1e8, 1e8;

  IDcontinuousPrior 				 initPrior(ipriors);
  lognormalPOLYNOMIAL_BMD_NC		 startingV(Y, X, suff_stat, 1);

  cBMDModel<lognormalPOLYNOMIAL_BMD_NC, IDcontinuousPrior> model(startingV, initPrior,
    isfixed1, fixed1, true);

  optimizationResult initoR = findMAP<lognormalPOLYNOMIAL_BMD_NC, IDcontinuousPrior>(&model);

  return initoR.max_parms;
}

// When analyzing lognormally distributed, summarized input data, transform
// arithmetic means and std. devs to their geometric equivalents.
// Note: This tranformation is approximate; the individual observations
// would be needed to calculate exact values.
void convertLognormalStats(vector <double>& mean, vector <double>& sd, int n)
{

  for (int i = 0; i < n; i++) {
    double x = mean[i];
    double y = sd[i];
    mean[i] = log(x) - log(1 + pow(y / x, 2.0)) / 2;
    sd[i] = sqrt(log(1.0 + pow(y / x, 2.0)));
  } // for

}

Eigen::MatrixXd getParmPriors(CModelID_t model, bool bLognormal, bool bConstVar,
                              int nparms, const Eigen::MatrixXd& max_parms,
                              int degree, bool isIncreasing, int restriction,
                              double yScale)
{
  Eigen::MatrixXd priors(nparms, NUM_PRIOR_COLS);

  if (bConstVar) {

    if (bLognormal) {
      // Lognormally distributed (implies constant variance)
      switch (model) {
        int j;
      case eExp2:
        expStashPriors[0] = max_parms(0, 0); // a
        expStashPriors[1] = fabs(log(fabs(max_parms(1,0))));           // b
        //expStashPriors[1] = bParm;           // b
        expStashPriors[2] = 7;               // c
        expStashPriors[3] = 2.0;              // d
        expStashPriors[4] = max_parms(2, 0); // log-alpha
        /* break; We want this to fall through */
      case eExp3: // Exp 3-5 are treated the same
      case eExp4:
      case eExp5:
        priors << 0, expStashPriors[0], 1, 0, 1e6,	//a
          0, expStashPriors[1], 1, 0, 18,	//b
          0, expStashPriors[2], 1, 0, 18,  //c
          0, expStashPriors[3], 1, 1, 18,  //d
          0, expStashPriors[4], 1, -18, 18; //log-alpha
        // Change c prior if downward direction of adversity
        if (!isIncreasing) {
          priors(2, 1) = -1; priors(2, 3) = -18;  priors(2, 4) = 0;
        }
        break;
      case eHill:
        priors <<
          0, fabs(max_parms(0, 0)), 1, 1e-8, 1e8, // g
          0, max_parms(1, 0), 1, -1e8, 1e8,      // v       
          0, 0.5, 1, 0, 100, // k
          0, 2, 0, 1e-8, 1000,  // n
          0, max_parms(2, 0), 1, -1000, 1000; // log-alpha
        if (restriction > 0) priors(3, 3) = 1;
          break;
      case ePoly:
        // beta0, beta1 priors
        priors(0, 0) = 0; priors(0, 1) = max_parms(0, 0); priors(0, 2) = 1; priors(0, 3) = -1e8; priors(0, 4) = 1e8;
        priors(1, 0) = 0; priors(1, 1) = max_parms(1, 0); priors(1, 2) = 1; priors(1, 3) = -1e8; priors(1, 4) = 1e8;
        // beta2 - betaN priors
        for (int i = 2; i <= degree; i++) {
          priors(i, 0) = 0; priors(i, 1) = 0; priors(i, 2) = 1; priors(i, 3) = -1e8; priors(i, 4) = 1e8;
        }
        // log-alpha prior
        j = nparms - 1;
        priors(j, 0) = 0; priors(j, 1) = max_parms(2, 0); priors(j, 2) = 1; priors(j, 3) = -1000; priors(j, 4) = 1000;
        if (restriction > 0) for (int i = 1; i <= degree; i++) priors(i, 3) = 0;
        else if (restriction < 0) for (int i = 1; i <= degree; i++) priors(i, 4) = 0;
        break;
      case ePow:
        priors <<
          0, max_parms(0, 0), 1, 1e-8, 1e8,    // g
          0, max_parms(1, 0), 1, -1e8, 1e8,    // b
          0, 1, 1, 1e-8, 100,                  // n
          0, max_parms(2, 0), 1, -1000, 1000;  // log-alpha
        if (restriction > 0) priors(2, 3) = 1;
        break;
      default:
        break;
      } // end switch

    }
    else {
      // Normally distributed with constant variance
      switch (model) {
        int j;
      case eExp2:
        cout << "\n)()()( max_parms(0, 0)= "<<max_parms(0, 0) 
          <<", log(fabs(max_parms(1,0)))= " << log(fabs(max_parms(1,0)))
          //<< ", bParm= "<< bParm<<")()()(\n"
          << endl;
        expStashPriors[0] = 1; // max_parms(0, 0); // a
		expStashPriors[1] = log(fabs(max_parms(1,0)));           // b
        //expStashPriors[1] = bParm;  // b
        //expStashPriors[2] = 7;      // c
        expStashPriors[2] = 4;      // c
        //expStashPriors[3] = 2;      // d
        expStashPriors[3] = 1;      // d
        //expStashPriors[4] = max_parms(2, 0); // alpha
        expStashPriors[4] = max_parms(2, 0); // -log(yScale*yScale); // alpha
        //expStashPriors[0] = max_parms(0, 0); // a
        //expStashPriors[1] = 1; // b
        //expStashPriors[2] = 4; // c
        //expStashPriors[3] = 1; // d
        //expStashPriors[4] = 0; // alpha
        /* break; We want this to fall through */
      case eExp3: // Exp 3-5 are treated the same
      case eExp4:
      case eExp5:
        priors << 0, expStashPriors[0], 1, 0, 1e6,	// a
          0, expStashPriors[1], 1, 0, 18,	// b
          0, expStashPriors[2], 1, 0, 18,  // log(c)
          0, expStashPriors[3], 1, 1, 18,  // d
          0, expStashPriors[4], 1, -18, 18; // log(alpha)
        if (!isIncreasing) {
          priors(2, 1) = -4; priors(2, 3) = -18;  priors(2, 4) = 0;
        }
        break;
      case eHill:
        priors <<
          0, max_parms(0, 0), 1, -1e8, 1e8, // g
          0, max_parms(1, 0), 1, -1e8, 1e8, // v
          0, 0.5, 1, 0, 30,  // k
          0, 1, 0, 1e-8, 18, // n
          0, 0, 1, -1000, 1000; // log-alpa
        if (restriction > 0) priors(3, 3) = 1;
        break;
      case ePoly:
        // beta0, beta1 priors
				priors(0, 0) = 0; priors(0, 1) = fabs(max_parms(0, 0)); priors(0, 2) = 1; priors(0, 3) = -1e8; priors(0, 4) = 1e8;
        priors(1, 0) = 0; priors(1, 1) = max_parms(1, 0); priors(1, 2) = 1; priors(1, 3) = -1e8; priors(1, 4) = 1e8;
        // beta2 - betaN priors
        for (int i = 2; i <= degree; i++) {
          priors(i, 0) = 0; priors(i, 1) = 0; priors(i, 2) = 1; priors(i, 3) = -1e8; priors(i, 4) = 1e8;
        }
        // variance coefficient prior
        j = nparms - 1;
        priors(j, 0) = 0; priors(j, 1) = max_parms(2, 0); priors(j, 2) = 1; priors(j, 3) = -1000; priors(j, 4) = 1000;
        if (restriction > 0) for (int i = 1; i <= degree; i++) priors(i, 3) = 0;
        else if (restriction < 0) for (int i = 1; i <= degree; i++) priors(i, 4) = 0;
        break;
      case ePow:
        priors <<
          0, max_parms(0, 0), 1, 1e-8, 1e8,
          0, max_parms(1, 0), 1, -1e8, 1e8,
          0, 1, 1, 1e-8, 100,
          0, max_parms(2, 0), 1, -1000, 1000;
        if (restriction > 0) priors(2, 3) = 1;

        break;
      default:
        break;
      } // end switch

    } // endif (bLognormal)
  }
  else {
    // Normally distributed with modeled variance
    switch (model) {
      int j;
    case eExp2:
      expStashPriors[0] = 1; // max_parms(0, 0); // a
	  expStashPriors[1] = log(fabs(max_parms(1,0))); // b
      //expStashPriors[1] = bParm;           // b
      expStashPriors[2] = 4; // c
      //expStashPriors[3] = 2; // d
      expStashPriors[3] = 1; // d
      expStashPriors[4] = 0; // rho
      //expStashPriors[5] = 0; // log-alpha
      expStashPriors[5] = max_parms(2, 0); // -log(yScale*yScale)
        //+ expStashPriors[4] * log(yScale); // alpha
    //expStashPriors[0] = max_parms(0, 0); // a
      //expStashPriors[1] = 1; // b
      //expStashPriors[2] = 4; // c
      //expStashPriors[3] = 2; // d
      //expStashPriors[4] = 0; // rho
      //expStashPriors[5] = 0; // log-alpha
     /* break; We want this to fall through */
    case eExp3: // Exp 3-5 are treated the same
    case eExp4:
    case eExp5:
      priors << 0, expStashPriors[0], 1, 0, 1e6, //a
        0, expStashPriors[1], 1, 0, 18, //b
        0, expStashPriors[2], 1, 0, 18, // log(c)
        0, expStashPriors[3], 1, 1, 18, //d
        0, expStashPriors[4], 1, -18, 18,//rho
        0, expStashPriors[5], 1, -18, 18; //log(alpha)
      if (!isIncreasing) {
        priors(2, 1) = -4; priors(2, 3) = -18;  priors(2, 4) = 0;
      }
      break;
    case eHill:
      priors << 0, max_parms(0, 0), 1, 1e-8, 1e8,
        0, max_parms(1, 0), 1, -1000, 1000,
        0, 0.5, 1, 0, 30,
        0, 1, 0, 1e-8, 18,
        0, 0, 1, -1000, 1000,
        0, 0, 1, -1000, 1000;
      if (restriction > 0) priors(3, 3) = 1;
      break;
    case ePoly:
      // beta0, beta1 priors
				  priors(0, 0) = 0; priors(0, 1) = fabs(max_parms(0, 0)); priors(0, 2) = 1; priors(0, 3) = -1e8; priors(0, 4) = 1e8;
      priors(1, 0) = 0; priors(1, 1) = max_parms(1, 0); priors(1, 2) = 1; priors(1, 3) = -1e8; priors(1, 4) = 1e8;
      // beta2 - betaN priors
      for (int i = 2; i <= degree; i++) {
        priors(i, 0) = 0; priors(i, 1) = 0; priors(i, 2) = 1; priors(i, 3) = -1e8; priors(i, 4) = 1e8;
      }
      // rho prior
      j = nparms - 2;
      priors(j, 0) = 0; priors(j, 1) = max_parms(2, 0); priors(j, 2) = 1; priors(j, 3) = -1000; priors(j, 4) = 1000;
      // variance coefficient prior
      j = nparms - 1;
      priors(j, 0) = 0; priors(j, 1) = max_parms(3, 0); priors(j, 2) = 1; priors(j, 3) = -1000; priors(j, 4) = 1000;
      if (restriction > 0) for (int i = 1; i <= degree; i++) priors(i, 3) = 0;
      else if (restriction < 0) for (int i = 1; i <= degree; i++) priors(i, 4) = 0;
      break;
    case ePow:
        priors << 0, fabs(max_parms(0, 0)), 1, 1e-8, 1e8,
        0, max_parms(1, 0), 1, -1e8, 1e8,
        0, 1, 1, 1e-8, 100,
        0, max_parms(2, 0), 1, -1000, 1000, // rho
        0, max_parms(3, 0), 1, -1000, 1000; // log-alpha
      if (restriction > 0) priors(2, 3) = 1;
      break;
    default:
      break;
    } // end switch

  } // endif (bConstVar)

  return priors;
}

#if 0
// Returns number of estimated parameters for the model
int modelDF(bmd_analysis *fitout, Eigen::MatrixXd priors,
            vector<bool>& bFixed, bool *zBounded) {
  int nparms = priors.rows();
  int estParmCount = nparms;
  for (int i = 0; i < nparms; i++) {
    zBounded[i] = false;
    double v = fitout->MAP_ESTIMATE(i, 0);
    if (bFixed[i]
        || fabs(v - priors(i, 3)) < BMDS_EPS
        || fabs(priors(i, 4) - v) < BMDS_EPS) {
      --estParmCount;
      zBounded[i] = true;
    } // end if
  } // end for
  return estParmCount;
} // modelDF
#endif // 0

// This routine gets the log-likelihoods for the data model tests and
// performs tests of interest (i.e., deviance) for the continuous models.
void analysisOfDeviance(ContinuousDeviance_t *zOut,
                        Eigen::MatrixXd& Y, Eigen::MatrixXd& X, double llFit,
                        int nEstParms, bool bLognormal, bool bConstVar,
                        bool bSuffStat, int dgCount) {

#ifndef R_COMPILATION
#if defined(_DEBUG) || defined(DEBUGLOG)
  ofstream file;
  file.open(LOGFILENAME);
#endif
#endif

  LLRow_t *llRow = NULL;
  // nObs = count of unique dose groups
  int nObs = dgCount; // X.rows();
  int nParms;
  int i; // loop variable

  //file << "nObs= " << nObs << " nEstParms= " << nEstParms << endl;


  // Calculate the additive constant value, which is included in BMDS 3.0
  // log-likelihoods but omitted by BMDS 2.x.
  //double additiveConstant = 0;
  //additiveConstant = (bSuffStat) ? -Y.col(2).sum() * log(2 * M_PI) / 2
  //                               : -X.rows() * log(2 * M_PI) / 2;

  // First, get the LLs for the data models
  if (bLognormal) {
    // *****  A1 Model  *****//
    llRow = &zOut->llRows[LK_A1];
    llRow->model = LK_A1;
    lognormalLLTESTA1 a1Test(Y, X, bSuffStat);
    nParms = a1Test.nParms();
    std::vector<double> fix1(nParms); for (i = 0; i < nParms; i++) { fix1[i] = 0.0; }
    std::vector<bool> isfix1(nParms); for (i = 0; i < nParms; i++) { isfix1[i] = false; }
    Eigen::MatrixXd a1Priors(nParms, NUM_PRIOR_COLS);
    for (i = 0; i < nParms; i++) {
      a1Priors.row(i) << 0, 0, 1, 0, 1e8;
    }
    a1Priors.row(nParms - 1) << 0, 0, 1, -1e8, 1e8;
    IDcontinuousPrior a1Init(a1Priors);
    statModel<lognormalLLTESTA1, IDcontinuousPrior> a1Model(a1Test, a1Init,
                                                            isfix1, fix1);
    optimizationResult a1Result = findMAP<lognormalLLTESTA1, IDcontinuousPrior>(&a1Model);
    llRow->ll = -a1Result.functionV; // +additiveConstant;
    //llRow->nParms = nObs + 1; // This is really the DF
    llRow->nParms = nParms; // This is really the DF
    llRow->aic = -2*llRow->ll + 2*llRow->nParms;
    //file << "A1 Model log-likelihood: " << llRow->ll << endl;

    // *****  A2 Model  *****//
    llRow = &zOut->llRows[LK_A2];
    llRow->model = LK_A2;
    lognormalLLTESTA2 a2Test(Y, X, bSuffStat);
    nParms = a2Test.nParms();
    std::vector<double> fix2(nParms); for (i = 0; i < nParms; i++) { fix2[i] = 0.0; }
    std::vector<bool> isfix2(nParms); for (i = 0; i < nParms; i++) { isfix2[i] = false; }
    Eigen::MatrixXd a2Priors(nParms, NUM_PRIOR_COLS);
    for (i = 0; i < nParms; i++) {
      a2Priors.row(i) << 0, 0, 1, -1e8, 1e8;
    }
    for (i = 0; i < nParms / 2; i++) {
      a2Priors.row(i) << 0, a1Result.max_parms(i), 1, 0, 1e8;
    }
    IDcontinuousPrior a2Init(a2Priors);
    statModel<lognormalLLTESTA2, IDcontinuousPrior> a2Model(a2Test, a2Init,
                                                            isfix2, fix2);
    optimizationResult a2Result = findMAP<lognormalLLTESTA2, IDcontinuousPrior>(&a2Model);
    llRow->ll = -a2Result.functionV; // +additiveConstant;
    //llRow->nParms = nObs * 2; // This is really the DF
    llRow->nParms = nParms; // This is really the DF
    llRow->aic = -2 * llRow->ll + 2 * llRow->nParms;
    //file << "A2 Model log-likelihood: " << llRow->ll << endl;

    // *****  A3 Model  *****//
    llRow = &zOut->llRows[LK_A3];
    llRow->model = LK_A3;
    // A3 = A1 for constant variance, which is always the case for lognormal
    llRow->ll = zOut->llRows[LK_A1].ll;
    llRow->nParms = zOut->llRows[LK_A1].nParms;
    llRow->aic = zOut->llRows[LK_A1].aic;
    //file << "A3 Model log-likelihood: " << llRow->ll << endl;

    // *****  Reduced Model  *****//
    llRow = &zOut->llRows[LK_R];
    llRow->model = LK_R;
    lognormalLLTESTR rTest(Y, X, bSuffStat);
    nParms = rTest.nParms();
    std::vector<double> fix4(nParms); for (i = 0; i < nParms; i++) { fix4[i] = 0.0; }
    std::vector<bool> isfix4(nParms); for (i = 0; i < nParms; i++) { isfix4[i] = false; }
    Eigen::MatrixXd rPriors(nParms, NUM_PRIOR_COLS);
    for (i = 0; i < nParms; i++) {
      rPriors.row(i) << 0, 1, 1, 0, 1e8;
    }
    rPriors.row(nParms - 1) << 0, 0, 1, -1e8, 1e8;
    IDcontinuousPrior rInit(rPriors);
    statModel<lognormalLLTESTR, IDcontinuousPrior> rModel(rTest, rInit,
                                                          isfix4, fix4);
    optimizationResult rResult = findMAP<lognormalLLTESTR, IDcontinuousPrior>(&rModel);
    llRow->ll = -rResult.functionV; // +additiveConstant;
    //llRow->nParms = 2; // This is really the DF
    llRow->nParms = nParms; // This is really the DF
    llRow->aic = -2 * llRow->ll + 2 * llRow->nParms;
    //file << "Reduced Model log-likelihood: " << llRow->ll << endl;

    // *****  Fitted Model  *****//
    llRow = &zOut->llRows[LK_FIT];
    zOut->llRows[LK_FIT].model = LK_FIT;
    zOut->llRows[LK_FIT].ll = llFit;
    zOut->llRows[LK_FIT].nParms = nEstParms;
    zOut->llRows[LK_FIT].aic = -2*llFit + 2*nEstParms;
    //file << "Fitted ll= " << llFit << " aic= " << zOut->llRows[LK_FIT].aic << endl;
  }
  else { // Normal distribution
    // *****  A1 Model  *****//
    llRow = &zOut->llRows[LK_A1];
    llRow->model = LK_A1;
    normalLLTESTA1 a1Test(Y, X, bSuffStat);

#ifdef DEBUG_HEAP
    cout << "***** Calling checkMemory @ line " << __LINE__ << endl;
    _CrtCheckMemory();
    //_CrtMemCheckpoint(&memSnapshot4); // Memory debugging
    //cout << "***** memDiff 3->4..." << endl;
    //if (_CrtMemDifference(&memStateDiff, &memSnapshot3, &memSnapshot4)) {
    //  _CrtMemDumpStatistics(&memStateDiff);
    //}
#endif // DEBUG_HEAP

    nParms = a1Test.nParms();
    std::vector<double> fix1(nParms); for (i = 0; i < nParms; i++) { fix1[i] = 0.0; }
    std::vector<bool> isfix1(nParms); for (i = 0; i < nParms; i++) { isfix1[i] = false; }
    Eigen::MatrixXd a1Priors(nParms, NUM_PRIOR_COLS);
    for (i = 0; i < nParms; i++) {
      a1Priors.row(i) << 0, 0, 1, -1e8, 1e8;
    }
    IDcontinuousPrior a1Init(a1Priors);
    statModel<normalLLTESTA1, IDcontinuousPrior> a1Model(a1Test, a1Init,
                                                            isfix1, fix1);
    optimizationResult a1Result = findMAP<normalLLTESTA1, IDcontinuousPrior>(&a1Model);
    llRow->ll = -a1Result.functionV; // +additiveConstant;
    //cout << "A1 nParms= " << nParms << ", LL= " << a1Result.functionV << endl;
    //llRow->nParms = nObs + 1; // This is really the DF
    llRow->nParms = nParms; // This is really the DF
    llRow->aic = -2 * llRow->ll + 2 * llRow->nParms;
    //file << "A1 Model log-likelihood: " << llRow->ll << endl;

    // *****  A2 Model  *****//
    llRow = &zOut->llRows[LK_A2];
    llRow->model = LK_A2;
    normalLLTESTA2 a2Test(Y, X, bSuffStat);
    nParms = a2Test.nParms();
    std::vector<double> fix2(nParms); for (i = 0; i < nParms; i++) { fix2[i] = 0.0; }
    std::vector<bool> isfix2(nParms); for (i = 0; i < nParms; i++) { isfix2[i] = false; }
    Eigen::MatrixXd a2Priors(nParms, NUM_PRIOR_COLS);
    for (i = 0; i < nParms; i++) {
      a2Priors.row(i) << 0, 0, 1, -1e8, 1e8;
    }
    for (i = 0; i < nParms / 2; i++) {
      a2Priors.row(i) << 0, a1Result.max_parms(i), 1, -1e8, 1e8;
    }
    IDcontinuousPrior a2Init(a2Priors);
    statModel<normalLLTESTA2, IDcontinuousPrior> a2Model(a2Test, a2Init,
                                                            isfix2, fix2);
    optimizationResult a2Result = findMAP<normalLLTESTA2, IDcontinuousPrior>(&a2Model);
    llRow->ll = -a2Result.functionV; // +additiveConstant;
    //cout << "A2 nParms= " << nParms << ", LL= " << a2Result.functionV << endl;
    //llRow->nParms = nObs * 2; // This is really the DF
    llRow->nParms = nParms; // This is really the DF
    llRow->aic = -2 * llRow->ll + 2 * llRow->nParms;
    //file << "A2 Model log-likelihood: " << llRow->ll << endl;

    // *****  A3 Model  *****//
    llRow = &zOut->llRows[LK_A3];
    llRow->model = LK_A3;
    if (bConstVar) { // A3 = A1 for constant variance
      llRow->ll = zOut->llRows[LK_A1].ll;
      llRow->nParms = zOut->llRows[LK_A1].nParms;
      llRow->aic = zOut->llRows[LK_A1].aic;
    }
    else {
      normalLLTESTA3 a3Test(Y, X, bSuffStat);
      nParms = a3Test.nParms();
      std::vector<double> fix3(nParms); for (i = 0; i < nParms; i++) { fix3[i] = 0.0; }
      std::vector<bool> isfix3(nParms); for (i = 0; i < nParms; i++) { isfix3[i] = false; }
      Eigen::MatrixXd a3Priors(nParms, NUM_PRIOR_COLS);
      for (i = 0; i < nParms; i++) {
        a3Priors.row(i) << 0, 0, 1, -1e8, 1e8;
      }
      for (i = 0; i < nParms - 2; i++) {
        a3Priors.row(i) << 0, a2Result.max_parms(i), 1, 0, 1e8;
      }
      IDcontinuousPrior a3Init(a3Priors);
      statModel<normalLLTESTA3, IDcontinuousPrior> a3Model(a3Test, a3Init,
                                                           isfix3, fix3);
      optimizationResult a3Result = findMAP<normalLLTESTA3, IDcontinuousPrior>(&a3Model);
      llRow->ll = -a3Result.functionV; // +additiveConstant;
      //cout << "A3 nParms= " << nParms << ", LL= " << a3Result.functionV << endl;
      //llRow->nParms = nObs * 2; // This is really the DF
      llRow->nParms = nParms; // This is really the DF
      llRow->aic = -2 * llRow->ll + 2 * llRow->nParms;
      //file << "A3 Model log-likelihood: " << llRow->ll << endl;
    } // endif bConstVar

    // *****  Reduced Model  *****//
    llRow = &zOut->llRows[LK_R];
    llRow->model = LK_R;
    normalLLTESTR rTest(Y, X, bSuffStat);
    nParms = rTest.nParms();
    std::vector<double> fix4(nParms); for (i = 0; i < nParms; i++) { fix4[i] = 0.0; }
    std::vector<bool> isfix4(nParms); for (i = 0; i < nParms; i++) { isfix4[i] = false; }
    Eigen::MatrixXd rPriors(nParms, NUM_PRIOR_COLS);
    for (i = 0; i < nParms; i++) {
      rPriors.row(i) << 0, 1, 1, 0, 1e8;
    }
    rPriors.row(nParms - 1) << 0, 0, 1, -1e8, 1e8;
    IDcontinuousPrior rInit(rPriors);
    statModel<normalLLTESTR, IDcontinuousPrior> rModel(rTest, rInit,
                                                       isfix4, fix4);
    optimizationResult rResult = findMAP<normalLLTESTR, IDcontinuousPrior>(&rModel);
    llRow->ll = -rResult.functionV; // +additiveConstant;
    //cout << " R nParms= " << nParms << ", LL= " << rResult.functionV << endl;
    //llRow->nParms = 2; // This is really the DF
    llRow->nParms = nParms; // This is really the DF
    llRow->aic = -2 * llRow->ll + 2 * llRow->nParms;
    //file << "Reduced Model log-likelihood: " << llRow->ll << endl;

    // *****  Fitted Model  *****//
    llRow = &zOut->llRows[LK_FIT];
    llRow->model = LK_FIT;
    llRow->ll = llFit;
    llRow->nParms = nEstParms;
    llRow->aic = -2*llFit + 2*nEstParms;
    //file << "Fitted ll= " << llFit << " aic= " << llRow->aic << endl;
  } // endif

  // Now, perform the tests of interest
  double dev; // deviance
  double pvalue;
  TestRow_t *testRow = NULL;
  int df;

  // Test #1 - A2 vs. Reduced - does mean and/or variance differ across dose groups
  testRow = &zOut->testRows[TI_1];
  testRow->testNumber = TI_1 + 1;
  testRow->deviance = dev = 2 * (zOut->llRows[LK_A2].ll - zOut->llRows[LK_R].ll);
  testRow->df = df = zOut->llRows[LK_A2].nParms - zOut->llRows[LK_R].nParms;
  testRow->pvalue = (dev < 0.0) ? -1 : 1.0 - gsl_cdf_chisq_P(dev, df);
#ifndef R_COMPILATION
#if defined(_DEBUG) || defined(DEBUGLOG)
  file << "test 1: llA2= " << zOut->llRows[LK_A2].ll << ", llR= " << zOut->llRows[LK_R].ll << ", dev= " << dev << ", df= " << df << ", pval= " << testRow->pvalue << endl;
#endif
#endif

  // Test #2 - A1 vs. A2 - homogeneity of variance across dose groups
  testRow = &zOut->testRows[TI_2];
  testRow->testNumber = TI_2 + 1;
  testRow->deviance = dev = 2 * (zOut->llRows[LK_A2].ll - zOut->llRows[LK_A1].ll);
  testRow->df = df = zOut->llRows[LK_A2].nParms - zOut->llRows[LK_A1].nParms;
  testRow->pvalue = (dev < 0.0) ? -1 : 1.0 - gsl_cdf_chisq_P(dev, df);
#ifndef R_COMPILATION
#if defined(_DEBUG) || defined(DEBUGLOG)
  file << "test 2: llA2= " << zOut->llRows[LK_A2].ll << ", llA1= " << zOut->llRows[LK_A1].ll << ", dev= " << dev << ", df= " << df << ", pval= " << testRow->pvalue << endl;
#endif
#endif

  // Test #3 - A2 vs. A3 - does the model describe the variances adequately
  testRow = &zOut->testRows[TI_3];
  testRow->testNumber = TI_3 + 1;
  testRow->deviance = dev = 2 * (zOut->llRows[LK_A2].ll - zOut->llRows[LK_A3].ll);
  testRow->df = df = zOut->llRows[LK_A2].nParms - zOut->llRows[LK_A3].nParms;
  testRow->pvalue = (dev < 0.0) ? -1 : 1.0 - gsl_cdf_chisq_P(dev, df);
#ifndef R_COMPILATION
#if defined(_DEBUG) || defined(DEBUGLOG)
  file << "test 3: llA2= " << zOut->llRows[LK_A2].ll << ", llA3 " << zOut->llRows[LK_A3].ll << ", dev= " << dev << ", df= " << df << ", pval= " << testRow->pvalue << endl;
#endif
#endif

  // Test #4 - A3 vs. Fitted model - does the fitted model describe the obs data adequately
  testRow = &zOut->testRows[TI_4];
  testRow->testNumber = TI_4 + 1;
  testRow->deviance = dev = 2 * (zOut->llRows[LK_A3].ll - llFit);
  testRow->df = df = zOut->llRows[LK_A3].nParms - nEstParms;
  testRow->pvalue = (dev < 0.0) ? -1 : 1.0 - gsl_cdf_chisq_P(dev, df);
#ifndef R_COMPILATION
#if defined(_DEBUG) || defined(DEBUGLOG)
  file << "test 4: llA3= " << zOut->llRows[LK_A3].ll << ", llFit= " << llFit << ", dev= " << dev << ", df= " << df << ", pval= " << testRow->pvalue << endl;
#endif
#endif

cleanup:
#ifndef R_COMPILATION
#if defined(_DEBUG) || defined(DEBUGLOG)
  file.close();
#endif
#endif

  return;
}

double calcInitSlopeValue(Eigen::MatrixXd X, Eigen::MatrixXd Y, bool bLognormal)
{
  double retValue = 0;
  cout << "calculating EXP model initial b parameter:" << endl;
  int nrows = X.rows();
  double xBar = X.col(0).mean(), yBar = Y.col(0).mean();
  cout << "xBar= " << xBar << ", yBar= " << yBar << endl;
  double ssx = 0, ssxy = 0;
  for (int i = 0; i < nrows; i++) {
    ssx += Y(i, 2) * (X(i, 0) - xBar) * (X(i, 0) - xBar);
    if (Y(i, 0) > 0) ssxy += Y(i, 2) * (X(i, 0) - xBar) * log(Y(i, 0));
    cout << "Loop " << i << ": ssx= " << ssx << ", ssxy= " << ssxy << endl;
  } // for
  retValue = fabs(ssxy / ssx);

  std::cout << "initial slope value= " << retValue << endl << endl; flush(std::cout);
  return(retValue);
}
