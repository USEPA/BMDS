// bmds_models.cpp : Defines the exported functions for the DLL application.
//
#ifndef _WINDOWS
#  include <unistd.h>
#endif // !WINDOWS
#include "stdafx.h" // Precompiled header - does nothing if building R version

#include <fstream>
#include <chrono>
#include <iostream>
#include <time.h>

#ifdef R_COMPILATION
    //necessary things to run in R    
    #include <RcppEigen.h>
    #include <RcppGSL.h>
#else 
	#include <Eigen/Dense>
#endif

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_randist.h>
    
#include <iomanip>


#include <iomanip>

#include <dBMDstatmod.h>
#include <log_likelihoods.h>
#include <binomModels.h>

#include "DichHillBMD_NC.h"
#include "DichMultistageBMD_NC.h"
#include "DichLogLogisticBMD_NC.h"
#include "DichLogProbitBMD_NC.h"
#include "DichWeibullBMD_NC.h"
#include "DichGammaBMD_NC.h"
#include "DichQlinearBMD_NC.h"
#include "DichLogisticBMD_NC.h"
#include "DichProbitBMD_NC.h"

#include "bmds_entry.h"
#include "bmd_calculate.h"
#include "binomialTests.h"
#include "bmds_dmodels.h"
#include "dModel.h"
//#include <math.h>

// Transform background parameters, which were computed using a logistic dist.
#define TRANSFORM_BG(x) 1.0 / (1.0 + exp(-(x)))
#define SET_FAILURE_RETURN(x) \
        (x)->MAP = 0;         \
        (x)->nparms = 0;      \
        (x)->BMD = -9999;     \
        (x)->BMDL = -9999;    \
        (x)->BMDU = -9999

#if defined WIN32 || defined _WINDOWS
#   ifndef R_COMPILATION
#       define ALLOC_MODEL_ID(x) SysAllocString(L ## x)
#   else
#       define ALLOC_MODEL_ID(x) x
#   endif      
#else
//TO LOUIS:
// DOESN'T THIS CODE CREATE A MEMORY LEAK?????????
//????????????????????????????????????????????????
//????????????????????????????????????????????????
//#define ALLOC_MODEL_ID(x) new char*(x) 
#define ALLOC_MODEL_ID(x) x
#endif // WIN32

// Starting value and step size for the analysis profile
constexpr auto ANALYSIS_ALPHA_VALUE = 0.005; // Starting value for generating profile
constexpr auto ANALYSIS_DEFAULT_STEP_SIZE = 0.03; // Step size for generating profile
constexpr auto MA_ANALYSIS_ALPHA_VALUE = 0.005; // MA starting value
constexpr auto MA_ANALYSIS_DEFAULT_STEP_SIZE = 0.03; // MA step size

// Starting value and step size to generate CDF values for BMD
constexpr auto BMD_CDF_START = 0.01;
constexpr auto BMD_CDF_STEP = 0.01;

using namespace std;

void _stdcall bmd_MA(BMDSInputType_t inputType,	 BMDSInputData_t *dataIn,
  MA_PRIORS *priors, double *p_m_probs,
  BMDS_D_Opts1_t *opt1, BMDS_D_Opts2_t *opt2,
  int n, double * post_p,
  double *ma_bmd, double *bmd, double *bmdl, double *bmdu) {

  gsl_set_error_handler_off();

  // Set default return values
  ma_bmd[0] = ma_bmd[1] = ma_bmd[2] = -9999;
  double background = opt1->background;;
  Eigen::MatrixXd  Y(n, 2); Eigen::MatrixXd  X(n, 1);
  for (int i = 0; i < n; i++) {
    Y(i, 0) = dataIn[i].response; X(i, 0) = dataIn[i].dose;
    Y(i, 1) = dataIn[i].groupSize;
  }
 
#ifdef DEBUGLOG
  ofstream file;
  file.open(LOGFILENAME);
  file << "bmrType= " << opt2->bmrType << " bmr= " << opt1->bmr
       << " alpha= " << opt1->alpha << ", background= " << background << endl;
  file << "Response & Group size:\n";
  file << Y << endl << endl;
  file << "Dose:\n";
  file << X << "\n" << endl;
#endif	
  std::vector<bool> fixedB_2(2); for (int i = 0; i < fixedB_2.size(); i++) { fixedB_2[i] = false; }
  std::vector<double> fixedV_2(2); for (int i = 0; i < fixedV_2.size(); i++) { fixedV_2[i] = 0.0; }
  std::vector<bool> fixedB_3(3); for (int i = 0; i < fixedB_3.size(); i++) { fixedB_3[i] = false; }
  std::vector<double> fixedV_3(3); for (int i = 0; i < fixedV_3.size(); i++) { fixedV_3[i] = 0.0; }
  std::vector<bool> fixedB_4(4); for (int i = 0; i < fixedB_4.size(); i++) { fixedB_4[i] = false; }
  std::vector<double> fixedV_4(4); for (int i = 0; i < fixedV_4.size(); i++) { fixedV_4[i] = 0.0; }

  //Log-Logit Setup
  Eigen::MatrixXd prior_loglogit(3, 5);
  prior_loglogit << 1.0, priors->mean_loglogit[0], priors->sd_loglogit[0], -20.0, 20.0,
    1.0, priors->mean_loglogit[1], priors->sd_loglogit[1], -40.0, 40.0,
    2.0, priors->mean_loglogit[2], priors->sd_loglogit[2], 1e-4, 20.0;
  IDbinomPrior loglogit_prior(prior_loglogit);

  //Gamma Setup
  Eigen::MatrixXd prior_gamma(3, 5);
  prior_gamma << 1.0, priors->mean_gamma[0], priors->sd_gamma[0], -18.0, 18.0,
    2.0, priors->mean_gamma[1], priors->sd_gamma[1], 0.2, 20.0,
    2.0, priors->mean_gamma[2], priors->sd_gamma[2], 1e-4, 100.0;
  IDbinomPrior gamma_prior(prior_gamma);

  //Logistic Setup
	Eigen::MatrixXd prior_logit(2, 5);
	prior_logit << 1.0, priors->mean_logistic[0], priors->sd_logistic[0],   -20,   20,
		   		   2.0, priors->mean_logistic[1], priors->sd_logistic[1], 1e-12, 100;
	IDbinomPrior logistic_prior(prior_logit);

	//Probit Setup
	Eigen::MatrixXd prior_probit(2, 5);
	prior_probit <<  1.0, priors->mean_probit[0], priors->sd_probit[0], -8, 8,
					2.0, priors->mean_probit[1], priors->sd_probit[1], 0, 40;
	IDbinomPrior probit_prior(prior_probit);

	//Quantal Linear Setup
	Eigen::MatrixXd prior_qlinear(2, 5);
	prior_qlinear << 1.0, priors->mean_qlinear[0], priors->sd_qlinear[0], -20  , 20,
				   2.0, priors->mean_qlinear[1], priors->sd_qlinear[1], 0, 100.0;
	IDbinomPrior qlinear_prior(prior_qlinear);

	//Log-Probit Setup
	Eigen::MatrixXd prior_logprobit(3, 5);
	prior_logprobit << 1.0, priors->mean_logprobit[0], priors->sd_logprobit[0], -20.0, 20.0,
					   1.0, priors->mean_logprobit[1], priors->sd_logprobit[1], -8, 8,
					   2.0, priors->mean_logprobit[2], priors->sd_logprobit[2], 1e-4, 40.0;
	IDbinomPrior logprobit_prior(prior_logprobit);

    //Weibull Setup
    Eigen::MatrixXd prior_weibull(3, 5);
    prior_weibull << 1.0, priors->mean_weibull[0], priors->sd_weibull[0], -20.0, 20.0,
      2.0, priors->mean_weibull[1], priors->sd_weibull[1], 1e-4, 18.0,
      2.0, priors->mean_weibull[2], priors->sd_weibull[2], 1e-4, 100.0;
    IDbinomPrior weibull_prior(prior_weibull);

    //Multistage Setup
  int mst_degree = 2; // Default polynomial degree
	Eigen::MatrixXd prior_mult2(3, 5);
	prior_mult2    << 1.0, priors->mean_mult2[0], priors->sd_mult2[0],-20.0, 20.0,
					  2.0, priors->mean_mult2[1], priors->sd_mult2[1], 1e-4, 100.0,
					  2.0, priors->mean_mult2[2], priors->sd_mult2[2], 1e-4, 1e6;
	IDbinomPrior mult2_prior(prior_mult2);

	//Hill Setup
	Eigen::MatrixXd prior_hill(4, 5);
	prior_hill << 1.0, priors->mean_hill[0], priors->sd_hill[0], -40, 40,
				  1.0, priors->mean_hill[1], priors->sd_hill[1], -40, 40,
				  1.0, priors->mean_hill[2], priors->sd_hill[2], -40, 40,
				  2.0, priors->mean_hill[3], priors->sd_hill[3], 1e-8, 40;
	IDbinomPrior    hill_prior(prior_hill);

#ifdef DEBUGLOG
  file << "logistic priors:\n" << prior_logit << endl;
  file << "probit priors:\n" << prior_probit << endl;
  file << "quantal linear priors:\n" << prior_qlinear << endl;
  file << "log-logistic priors:\n" << prior_loglogit << endl;
  file << "log-probit priors:\n" << prior_logprobit << endl;
  file << "multistage priors:\n" << prior_mult2 << endl;
  file << "weibull priors:\n" << prior_weibull << endl;
  file << "gamma priors:\n" << prior_gamma << endl;
  file << "hill priors:\n" << prior_hill << endl << endl;
#endif // DEBUGLOG
  std::vector<double> prior_prob(9); for (int i = 0; i < prior_prob.size(); i++) { prior_prob[i] = p_m_probs[i]; }
  std::vector<double> post_prob(9);
  std::list<bmd_analysis> analyses;
  
  try {
#ifdef DEBUGLOG
    file << "Running model: Logistic... "; file.flush();
#endif
    analyses.push_back(bmd_analysis_DNC<dich_logisticModelNC, IDbinomPrior>(Y, X, prior_logit,
      fixedB_2, fixedV_2, 0, opt1->bmr, opt2->bmrType == eExtraRisk, MA_ANALYSIS_ALPHA_VALUE, MA_ANALYSIS_DEFAULT_STEP_SIZE));
#ifdef DEBUGLOG
    file << "done." << endl;
    file << "Running model: Probit... "; file.flush();
#endif
    analyses.push_back(bmd_analysis_DNC<dich_probitModelNC, IDbinomPrior>(Y, X, prior_probit,
      fixedB_2, fixedV_2, 0, opt1->bmr, opt2->bmrType == eExtraRisk, MA_ANALYSIS_ALPHA_VALUE, MA_ANALYSIS_DEFAULT_STEP_SIZE));
#ifdef DEBUGLOG
    file << "done." << endl;
    file << "Running model: Quantal Linear... "; file.flush();
#endif
    analyses.push_back(bmd_analysis_DNC<dich_qlinearModelNC, IDbinomPrior>(Y, X, prior_qlinear,
      fixedB_2, fixedV_2, 0, opt1->bmr, opt2->bmrType == eExtraRisk, MA_ANALYSIS_ALPHA_VALUE, MA_ANALYSIS_DEFAULT_STEP_SIZE));
#ifdef DEBUGLOG
    file << "done." << endl;
    file << "Running model: Log-logistic... "; file.flush();
#endif
    analyses.push_back(bmd_analysis_DNC<dich_loglogisticModelNC, IDbinomPrior>(Y, X, prior_loglogit,
      fixedB_3, fixedV_3, 0, opt1->bmr, opt2->bmrType == eExtraRisk, MA_ANALYSIS_ALPHA_VALUE, MA_ANALYSIS_DEFAULT_STEP_SIZE));
#ifdef DEBUGLOG
    file << "done." << endl;
    file << "Running model: Log-probit... "; file.flush();
#endif
    analyses.push_back(bmd_analysis_DNC<dich_logProbitModelNC, IDbinomPrior>(Y, X, prior_logprobit,
      fixedB_3, fixedV_3, 0, opt1->bmr, opt2->bmrType == eExtraRisk, MA_ANALYSIS_ALPHA_VALUE, MA_ANALYSIS_DEFAULT_STEP_SIZE));
#ifdef DEBUGLOG
    file << "done." << endl;
    file << "Running model: Multistage... "; file.flush();
#endif
    analyses.push_back(bmd_analysis_DNC<dich_multistageNC, IDbinomPrior>(Y, X, prior_mult2,
      fixedB_3, fixedV_3, mst_degree, opt1->bmr, opt2->bmrType == eExtraRisk, MA_ANALYSIS_ALPHA_VALUE, MA_ANALYSIS_DEFAULT_STEP_SIZE));
#ifdef DEBUGLOG
    file << "done." << endl;
    file << "Running model: Weibull... "; file.flush();
#endif
    analyses.push_back(bmd_analysis_DNC<dich_weibullModelNC, IDbinomPrior>(Y, X, prior_weibull,
      fixedB_3, fixedV_3, 0, opt1->bmr, opt2->bmrType == eExtraRisk, MA_ANALYSIS_ALPHA_VALUE, MA_ANALYSIS_DEFAULT_STEP_SIZE));
#ifdef DEBUGLOG
    file << "done." << endl;
    file << "Running model: Gamma... "; file.flush();
#endif
    analyses.push_back(bmd_analysis_DNC<dich_gammaModelNC, IDbinomPrior>(Y, X, prior_gamma,
      fixedB_3, fixedV_3, 0, opt1->bmr, opt2->bmrType == eExtraRisk, MA_ANALYSIS_ALPHA_VALUE, MA_ANALYSIS_DEFAULT_STEP_SIZE));
#ifdef DEBUGLOG
    file << "done." << endl;
    file << "Running model: Dichotomous Hill... "; file.flush();
#endif
    analyses.push_back(bmd_analysis_DNC<dich_hillModelNC, IDbinomPrior>(Y, X, prior_hill,
      fixedB_4, fixedV_4, 0, opt1->bmr, opt2->bmrType == eExtraRisk, MA_ANALYSIS_ALPHA_VALUE, MA_ANALYSIS_DEFAULT_STEP_SIZE));
#ifdef DEBUGLOG
    file << "done." << endl;
#endif
  } // try
  catch (std::exception &e) {
#ifdef DEBUGLOG
    file << "Error - analysis failed with an unexpected exception:" << endl;
    file << e.what();
    file.close();
#endif
    return;
  } // catch
  catch (...) {
#ifdef DEBUGLOG
    file << "Error - analysis failed with an unknown exception:" << endl;
    file.close();
#endif
    return;
  } // catch

	////////////////////////////////////////////////////////////////////////////
	// compute the log Laplace Approximation
	int j = 0;
	double v; 
	
	for (bmd_analysis b : analyses) {
		// note if the estimate of the covariance matrix is near singular, the estimate
		// for the determinant may be negative.  In this case we let it be zero
		post_prob[j] = b.MAP_ESTIMATE.rows() / 2 * log(2 * M_PI) - b.MAP + 0.5*log(max(0.0,b.COV.determinant()));
		v = post_prob[j]; 
		j++;
	}

	double max_prob = *max_element(post_prob.begin(), post_prob.end());
	double norm_sum = 0.0;
	
	for (int i = 0; i < post_prob.size(); i++) {
		post_prob[i] = post_prob[i] - max_prob + log(prior_prob[i]);
		norm_sum += exp(post_prob[i]);
	}

	for (int i = 0; i < post_prob.size(); i++) {
		post_prob[i] = exp(post_prob[i]) / norm_sum;
		post_p[i] = double(post_prob[i]);
		v = post_p[i]; 
	}

	//returnV->BMD = find_maBMDquantile(0.50, post_prob, analyses);
	float temp;
	temp = float(find_maBMDquantile(opt1->alpha, post_prob, analyses));
	ma_bmd[1] = double(temp); // find_maBMDquantile(opt1->alpha, post_prob, analyses);
	temp = float(find_maBMDquantile(1.0 - opt1->alpha, post_prob, analyses));
	ma_bmd[2] = double(temp); // find_maBMDquantile(1.0 - opt1->alpha, post_prob, analyses);
	int k = 0; 
	ma_bmd[0] = 0.0; 
	// individual BMD analysis
#ifdef DEBUGLOG
  file << "Setting result arrays for the individual models before exiting" << endl;
  file.flush();
#endif
	for (bmd_analysis n : analyses) {
		temp = (float(n.MAP_BMD)); 
		bmd[k] = double(temp);
		temp = float(post_prob[k] * n.MAP_BMD);
		ma_bmd[0] += double(temp);
		temp = float(n.BMD_CDF.inv(opt1->alpha));
		bmdl[k] = double(temp);
		temp = float(n.BMD_CDF.inv(1.0 - opt1->alpha));
		bmdu[k] = double(temp); // n.BMD_CDF.inv(1.0 - opt1->alpha);
		k++;
	}
	k = 0; 

  cleanup:
#ifdef DEBUGLOG
    file.close();
#endif
    return;
}

// This is the main entry point for running the models
int _stdcall run_dmodel2(DModelID_t *p_m_id, BMD_ANAL *returnV,
                        BMDSInputType_t *p_inputType, BMDSInputData_t *dataIn,
                        PRIOR *priorsIn, BMDS_D_Opts1_t *opt1, BMDS_D_Opts2_t *opt2, int *p_n) {
  int exitStatus = -1; // Default to error status
  DModelID_t m_id = *p_m_id;
  BMDSInputType_t inputType = *p_inputType;
  int n = *p_n;
  string myName(MODELNAME[m_id]);

  // Uncomment for Matt's debugging
  //ofstream file;
  //file.open("bob.out");

#ifdef DEBUGLOG
  ofstream file;
  string fname(__FUNCTION__);
  fname.append('.' + myName + ".log");
  file.open(fname, fstream::app);
  time_t now = time(0);
  struct tm localtm;
  char str[64];
  localtime_s(&localtm, &now);
  asctime_s(str, sizeof(str), &localtm);
  file << "Running model " << myName << " at " << str << endl;
  //file << "Function parameter addresses:" << endl;
  //file << std::hex;
  //file << "m_id= 0x" << p_m_id << " inputType= 0x" << p_inputType << endl;
  //file << "dataIn= 0x" << dataIn << " priorsIn= 0x" << priorsIn << " n= 0x" << p_n << endl;
  //file << "returnV= 0x" << returnV << " opt1= 0x" << opt1 << " opt2= 0x" << opt2 << endl;
  //file << std::dec;
  //file << "sizeof returnV= " << sizeof(returnV) << " sizeof opt1= " << sizeof(opt1)
  //  << " sizeof opt2= " << sizeof(opt2) << '\n' << endl;
  file.flush();
#endif // DEBUGLOG

  gsl_set_error_handler_off();
  dModel *pModel = NULL;
  int nparms = 0;
  int degree = 0; // Only used for multistage
  // Initialize return values
  returnV->model_id = ALLOC_MODEL_ID("1.1");
  returnV->MAP = 0;
  returnV->BMD = -9999;
  returnV->BMDL = -9999;
  returnV->BMDU = -9999;
  returnV->nparms = 0;
  bmd_analysis a;
  double background = opt1->background;
  double cl_alpha = opt1->alpha; // Confidence limit value from user
  double a_alpha = ANALYSIS_ALPHA_VALUE; // Starting value for profile
  double myVal; // Used to generate CDF values
  // Find the maximum dose to scale doses by the max
  double dMax = 0;
  for (int i = 0; i < n; i++) {
    if (dataIn[i].dose > dMax) dMax = dataIn[i].dose;
  }
  Eigen::MatrixXd  Y(n, 2); Eigen::MatrixXd  X(n, 1);
  for (int i = 0; i < n; i++) {
    Y(i, 0) = dataIn[i].response; X(i, 0) = dataIn[i].dose / dMax;
    Y(i, 1) = dataIn[i].groupSize;
  }
#ifdef DEBUGLOG
  file << "bmrType= " << opt2->bmrType << " bmr= " << opt1->bmr
    << " alpha= " << cl_alpha << ", background= " << background
    << " inputType= " << inputType << " n= " << n << '\n' << endl;
  file << "Input data with doses scaled by: " << dMax << endl;
  file << " Dose\tResp\tN" << endl;

  for (int i = 0; i < n; i++) {
    file << " " << X(i, 0)
      << '\t' << Y(i, 0)
      << '\t' << Y(i, 1)
      << '\n';
  } // for
  file << endl;
#endif

  switch (m_id) {
  case eDHill:
    nparms = NPARM_DHILL;
    pModel = new bmds_dhill();
    break;
  case eGamma:
    nparms = NPARM_GAMMA;
    pModel = new bmds_gamma();
    break;
  case eLogistic:
    nparms = NPARM_LOGISTIC;
    pModel = new bmds_logistic();
    break;
  case eLogLogistic:
    nparms = NPARM_LNLOGISTIC;
    pModel = new bmds_loglogistic();
    break;
  case eLogProbit:
    nparms = NPARM_LNPROBIT;
    pModel = new bmds_logprobit();
    break;
  case eMultistage:
    degree = opt2->degree;
    nparms = NPARM_MSTAGE(degree);
    myName.append('-' + to_string(degree));
    pModel = new bmds_multistage(degree);
    break;
  case eProbit:
    nparms = NPARM_PROBIT;
    pModel = new bmds_probit();
    break;
  case eQLinear:
    nparms = NPARM_QLINEAR;
    pModel = new bmds_qlinear();
    break;
  case eWeibull:
    nparms = NPARM_WEIBULL;
    pModel = new bmds_weibull();
    break;
  default:
    break;
  } // end switch
  Eigen::MatrixXd priors(nparms, NUM_PRIOR_COLS);
  std::vector<bool> fixedB(nparms);
  std::vector<double> fixedV(nparms);
  for (int i = 0; i < nparms; i++) {
    priors(i, 0) = priorsIn[i].type;
    priors(i, 1) = priorsIn[i].initalValue;
    priors(i, 2) = priorsIn[i].stdDev;
    priors(i, 3) = priorsIn[i].minValue;
    priors(i, 4) = priorsIn[i].maxValue;
    // Initialize fixed parameter vectors.
    fixedB[i] = false;
    fixedV[i] = 0.0;
  }
#ifdef DEBUGLOG
  file << "Parameter priors:\n" << priors << endl;
#endif

  if (background != BMDS_BLANK_VALUE) {
    if (m_id != eLogistic && m_id != eProbit) background = fixedParmToLogistic(background);
#ifdef DEBUGLOG
    file << "Fixing parameter# " << iBACKGROUND << " to " << background << endl;
#endif
    priors(0, 1) = fixedV[iBACKGROUND] = background;
    fixedB[iBACKGROUND] = true;
  }
  try {
    switch (m_id) {
    case eDHill:
      a = bmd_analysis_DNC<dich_hillModelNC, IDbinomPrior>(Y, X, priors,
                                                           fixedB, fixedV, degree,
                                                           opt1->bmr, opt2->bmrType == eExtraRisk,
                                                           a_alpha, ANALYSIS_DEFAULT_STEP_SIZE);
      break;
    case eGamma:
      a = bmd_analysis_DNC<dich_gammaModelNC, IDbinomPrior>(Y, X, priors,
                                                            fixedB, fixedV, degree,
                                                            opt1->bmr, opt2->bmrType == eExtraRisk,
                                                            a_alpha, ANALYSIS_DEFAULT_STEP_SIZE);
      break;
    case eLogistic:
      a = bmd_analysis_DNC<dich_logisticModelNC, IDbinomPrior>(Y, X, priors,
                                                               fixedB, fixedV, degree,
                                                               opt1->bmr, opt2->bmrType == eExtraRisk,
                                                               a_alpha, ANALYSIS_DEFAULT_STEP_SIZE);
      break;
    case eLogLogistic:
      a = bmd_analysis_DNC<dich_loglogisticModelNC, IDbinomPrior>(Y, X, priors,
                                                                  fixedB, fixedV, degree,
                                                                  opt1->bmr, opt2->bmrType == eExtraRisk,
                                                                  a_alpha, ANALYSIS_DEFAULT_STEP_SIZE);
       break;
    case eLogProbit:
      a = bmd_analysis_DNC<dich_logProbitModelNC, IDbinomPrior>(Y, X, priors,
                                                                fixedB, fixedV, degree,
                                                                opt1->bmr, opt2->bmrType == eExtraRisk,
                                                                a_alpha, ANALYSIS_DEFAULT_STEP_SIZE);
      break;
    case eMultistage:
      a = bmd_analysis_DNC<dich_multistageNC, IDbinomPrior>(Y, X,
                                                            priors, fixedB, fixedV, degree, opt1->bmr,
                                                            opt2->bmrType == eExtraRisk, a_alpha, ANALYSIS_DEFAULT_STEP_SIZE);
      break;
    case eProbit:
      a = bmd_analysis_DNC<dich_probitModelNC, IDbinomPrior>(Y, X, priors,
                                                             fixedB, fixedV, degree,
                                                             opt1->bmr, opt2->bmrType == eExtraRisk,
                                                             a_alpha, ANALYSIS_DEFAULT_STEP_SIZE);
      break;
    case eQLinear:
      a = bmd_analysis_DNC<dich_qlinearModelNC, IDbinomPrior>(Y, X, priors,
                                                              fixedB, fixedV, degree,
                                                              opt1->bmr, opt2->bmrType == eExtraRisk,
                                                              a_alpha, ANALYSIS_DEFAULT_STEP_SIZE);
      break;
    case eWeibull:
      a = bmd_analysis_DNC<dich_weibullModelNC, IDbinomPrior>(Y, X, priors,
                                                              fixedB, fixedV, degree,
                                                              opt1->bmr, opt2->bmrType == eExtraRisk,
                                                              a_alpha, ANALYSIS_DEFAULT_STEP_SIZE);
      break;
    default:
      break;
    } // end switch
  } // try
  catch (std::exception &e) {
#ifdef DEBUGLOG
    file << "Error - analysis failed with an unexpected exception:" << endl;
    file << e.what();
#endif
    goto cleanup;
  } // catch
  catch (...) {
#ifdef DEBUGLOG 
    file << "Error - analysis failed with an unknown exception:" << endl;
#endif
    goto cleanup;
  } // catch

#ifdef DEBUGLOG
  file << "Model Results BEFORE re-scaling:" << endl;
  file << std::scientific << setprecision(17);
  file << " MAP (-LL) = " << a.MAP << '\n';
  file << " BMD= " << a.MAP_BMD << " BMDL= " << a.BMD_CDF.inv(cl_alpha) << " BMDU= "
    << a.BMD_CDF.inv(1.0 - cl_alpha) << '\n';
  file << "MLE Parameter values\n";
  for (int i = 0; i < nparms; i++) {
    file << "Parm" << i << '\t' << a.MAP_ESTIMATE(i, 0) << endl;
  }
  file << "\n" << endl << flush;
#endif // DEBUGLOG

  // The BIC-equivalent is only valid for Bayesian runs, although we always
  // set it for now.
 
  for (int i = 0; i < a.COV.rows(); i++){
    for (int j = 0; j < a.COV.rows(); j++){
      returnV->covM[i + j*a.COV.rows()] = a.COV(i,j); 
    
    }
  }
  
  returnV->BIC_Equiv = a.MAP_ESTIMATE.rows() / 2 * log(2 * M_PI) - a.MAP + 0.5*log(max(0.0, a.COV.determinant()));
  returnV->MAP = -a.MAP; // a.MAP is the negative log-likelihood
  returnV->BMD = a.MAP_BMD * dMax;
  returnV->BMDL = a.BMD_CDF.inv(cl_alpha) * dMax;
  returnV->BMDU = a.BMD_CDF.inv(1.0 - cl_alpha) * dMax;
  returnV->nparms = nparms;

  // Populate CDF values to return
  myVal = BMD_CDF_START;
  for (int i = 0; i < returnV->nCDF; i++) {
    returnV->aCDF[i] = a.BMD_CDF.inv(myVal) * dMax;
    myVal += BMD_CDF_STEP;
  }
  // Rescale parameters affected by dose scaling
  // The background is treated the same across models. The other MLE values are
  // copied to the return struct, since many do not require rescaling. Those that
  // do are adjusted afterward.
  returnV->PARMS[0] = TRANSFORM_BG(a.MAP_ESTIMATE(0, 0));
  if (background != BMDS_BLANK_VALUE && returnV->PARMS[0] <= 1e-6) returnV->PARMS[0] = 0;
  for (int i = 1; i < nparms; i++) {
    returnV->PARMS[i] = a.MAP_ESTIMATE(i, 0);
  }
  
  // Correct/re-scale MLE values for models that need it
  switch (m_id) {
  case eDHill:
    // The N parameter was determined with a distribution similar to G's
    returnV->PARMS[1] = TRANSFORM_BG(a.MAP_ESTIMATE(1, 0));
    // The A parameter is affected by dose-scaling
    returnV->PARMS[2] = -log(dMax) * a.MAP_ESTIMATE(3, 0) + a.MAP_ESTIMATE(2, 0);
    break;
  case eGamma:
    returnV->PARMS[2] = a.MAP_ESTIMATE(2, 0) / dMax;
    break;
  case eLogLogistic:
  case eLogProbit:
    // The A parameter is affected by dose-scaling
    returnV->PARMS[1] = -log(dMax) * a.MAP_ESTIMATE(2, 0) + a.MAP_ESTIMATE(1, 0);
    break;
  case eLogistic:
  case eProbit:
    // background is NOT transformed for logistic or probit
    returnV->PARMS[0] = a.MAP_ESTIMATE(0, 0);
    // no break - we want to fall through
  case eQLinear:
    // Slope parameter is affected by dose-scaling
    returnV->PARMS[1] = a.MAP_ESTIMATE(1, 0) / dMax;
    break;
  case eMultistage:
    for (int i = 1; i < nparms; i++) {
      returnV->PARMS[i] = a.MAP_ESTIMATE(i, 0) / pow(dMax, i);
    }
    
    break;
  case eWeibull:
    returnV->PARMS[2] = a.MAP_ESTIMATE(2, 0) * pow(dMax, -a.MAP_ESTIMATE(1, 0));
    break;
  default:
    break;
  } // end switch

  int estParmCount; // Number of model parameters that are estimated
#ifdef DEBUGLOG
  file << "Model Results BEFORE re-scaling (AGAIN):" << endl;
  file << std::scientific << setprecision(17);
  file << " MAP (-LL) = " << a.MAP << '\n';
  file << " BMD= " << a.MAP_BMD << " BMDL= " << a.BMD_CDF.inv(cl_alpha) << " BMDU= "
    << a.BMD_CDF.inv(1.0 - cl_alpha) << '\n';
  file << "MLE Parameter values\n";
  for (int i = 0; i < nparms; i++) {
    file << "Parm" << i << '\t' << a.MAP_ESTIMATE(i, 0) << endl;
  }
  file << "\n" << endl << flush;
#endif // DEBUGLOG
  estParmCount = modelDF(&a, priors, fixedB, returnV->boundedParms);
  pModel->setParms(returnV->PARMS);
  pModel->calcGoF(returnV->gof, dataIn, n, estParmCount);
  // Calculate AIC value for frequentist runs
  double aic;
  // Multiply by +2 not -2 b/c a.MAP is already a positive value
  aic = 2 * a.MAP + 2 * estParmCount;
#ifdef DEBUGLOG
  file << "Estimated parameter count = " << estParmCount << ", AIC = " << aic << endl;
#endif
  returnV->AIC = aic;

  analysisOfDeviance(returnV->deviance, Y, X, -a.MAP, estParmCount);

#ifdef DEBUGLOG
  file << "Result values:" << endl;
  file << "LL= " << returnV->MAP << " bmd= " << returnV->BMD << " bmdl= "
    << returnV->BMDL << " bmdu= " << returnV->BMDU << endl;
  for (int i = 0; i < returnV->nparms; i++) {
    file << "  Parm" << i << '\t' << returnV->PARMS[i] << endl;
  }
  file << "\n" << endl << flush;
#endif
  exitStatus = 0; // Success

cleanup:
#ifdef DEBUGLOG
  file.close();
#endif
  
  return exitStatus;
} // end run_dmodel()

// Returns number of estimated parameters for the model.
//
int modelDF(bmd_analysis *fitout, Eigen::MatrixXd priors,
            vector<bool>& bFixed, bool *zBounded) {
  int nparms = priors.rows();
  int estParmCount = nparms;
#ifdef DEBUGLOG
  ofstream file;
  file.open(LOGFILENAME, fstream::app);
  char str[64]; memset(str, '\0', 64);
#ifdef _WINDOWS
  time_t now = time(0);
  struct tm localtm;
  localtime_s(&localtm, &now);
  asctime_s(str, sizeof(str), &localtm);
#endif // _WINDOWS
  file << "In modelDF "<< str << endl << flush;
  file << std::scientific << setprecision(17);
#endif
  for (int i = 0; i < nparms; i++) {
    zBounded[i] = false;
    double v = fitout->MAP_ESTIMATE(i, 0);
    DEBUG_LOG(file, "i=" << i << "\tparm= " << v << ", min= " << priors(i, 3) << ", max= " << priors(i, 4));
    if (bFixed[i]
        || fabs(v - priors(i, 3)) < BMDS_EPS
        || fabs(priors(i, 4) - v) < BMDS_EPS) {
      --estParmCount;
      zBounded[i] = true;
      DEBUG_LOG(file, "  * bounded = true");
    } // end if
  } // end for
  DEBUG_LOG(file, " ");
  DEBUG_CLOSE_LOG(file);
  return estParmCount;
} // modelDF

void analysisOfDeviance(DichotomousDeviance_t *zOut, Eigen::MatrixXd& Y,
                       Eigen::MatrixXd& X, double llFit, int nEstParms) {
  gsl_set_error_handler_off();

  // Initialize return values
  zOut->llFull = 0;
  zOut->llReduced = 0;
  zOut->devReduced = 0;
  zOut->pvReduced = 0;
#ifdef DEBUGLOG
  ofstream file;
  file.open(LOGFILENAME);
  file << "Dose\tResponse\tSize\n";
  file << X << Y << endl << flush;
#endif

  int nObs = X.rows();
  int df = nObs - nEstParms;
  int nParms = 0; // #Parms for current test model (full, fitted, reduced)
#ifdef DEBUGLOG
  file << "nObs= " << nObs << " nEstParms= " << nEstParms << endl << flush;
#endif
// Data structures for the full model
  binomialLLTESTA1 BTest1(Y, X);
  nParms = BTest1.nParms(); // Should == nObs
#ifdef DEBUGLOG
   file << "Full model #Parms = " << nParms << endl << flush;
#endif
  std::vector<double> fix1(nParms); for (int i = 0; i < fix1.size(); i++) { fix1[i] = 0.0; }
  std::vector<bool> isfix1(nParms); for (int i = 0; i < isfix1.size(); i++) { isfix1[i] = false; }
#ifdef DEBUGLOG  
	file << NUM_PRIOR_COLS << endl << flush; 
#endif
   Eigen::MatrixXd priors1(nParms, NUM_PRIOR_COLS);
  for (int i = 0; i < priors1.rows(); i++) {
    priors1.row(i) << 0, 0, 1, -20, 20;
  }
#ifdef DEBUGLOG
  file << priors1 << endl << flush; 
#endif 
  IDbinomPrior model_prior(priors1);
  statModel<binomialLLTESTA1, IDbinomPrior> saturated_model(BTest1, model_prior, isfix1, fix1);
  optimizationResult sat_mod_result = findMAP<binomialLLTESTA1, IDbinomPrior>(&saturated_model);
  zOut->llFull = -sat_mod_result.functionV;
  zOut->nparmFull = nParms;
#ifdef DEBUGLOG  
	file << "Saturated Model log-likelihood: " << zOut->llFull << endl << flush;
#endif  
  zOut->devFit = 2 * (zOut->llFull - llFit);
  zOut->dfFit = df;
  zOut->nparmFit = nEstParms;
#ifdef DEBUGLOG
  file << zOut->devFit << endl; 
#endif  
  zOut->pvFit = 1.0 - gsl_cdf_chisq_P(zOut->devFit, df);
#ifdef DEBUGLOG
  file << "Fitted model ll= " << llFit << " dev= " << zOut->devFit
       << " df= " << df << ", pv= " << zOut->pvFit << endl << flush;
#endif 
  // Data structures for the reduced model
  binomialLLTESTA2 BTest2(Y, X);
  nParms = BTest2.nParms(); // Should == 1
  std::vector<double> fix2(nParms); for (int i = 0; i < fix2.size(); i++) { fix2[i] = 0.0; }
  std::vector<bool> isfix2(nParms); for (int i = 0; i < isfix2.size(); i++) { isfix2[i] = false; }
  Eigen::MatrixXd priors2(nParms, NUM_PRIOR_COLS);
  for (int i = 0; i < priors2.rows(); i++) {
    priors2.row(i) << 0, 0, 1, -20, 20;
  }
  IDbinomPrior model_prior2(priors2);
  statModel<binomialLLTESTA2, IDbinomPrior> reduced_model(BTest2, model_prior2, isfix2, fix2);
  optimizationResult red_mod_result = findMAP<binomialLLTESTA2, IDbinomPrior>(&reduced_model);
  zOut->llReduced = -red_mod_result.functionV;
  zOut->devReduced = 2 * (zOut->llFull - zOut->llReduced);
  df = nObs - 1;
  zOut->dfReduced = df;
  zOut->nparmReduced = nParms;
  zOut->pvReduced = 1.0 - gsl_cdf_chisq_P(zOut->devReduced, df);
#ifdef DEBUGLOG
  file << "Reduced model ll= " << zOut->llReduced << " dev= " << zOut->devReduced
    << " df= " << df << ", pv= " << zOut->pvReduced << endl << flush;
#endif

cleanup:

#ifdef DEBUGLOG
  file.close();
#endif

  return;
}

double fixedParmToLogistic(double x) {
  double val;
  if (x <= 0) val = -18;
  else if (x >= 1) val = 18;
  else val = -log(1 / x - 1);
  return val;
}
