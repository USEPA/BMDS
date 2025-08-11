#ifdef __cplusplus
#  include <Eigen/Dense>
#endif
// #include "bmdStruct.h"
#include "dichotomous_entry_code.h"

// define __stdcall to nothing if not on Windows platform
#ifndef WIN32
#  define __stdcall
#endif

// define import/export attributes if on Windows platform
#ifndef R_COMPILATION
#  ifndef _WIN32
#    define BMDS_ENTRY_API
#  else
#    ifdef RBMDS_EXPORTS
#      define BMDS_ENTRY_API __declspec(dllexport)
#    else
#      define BMDS_ENTRY_API __declspec(dllimport)
#    endif  // BMDS_MODELS_EXPORTS
#  endif
#else
#  define BMDS_ENTRY_API
#endif

const double BMDS_EPS = 1.0e-6;
const double CloseToZero = 1.0e-7;
const double Log_zero = -22;  // e^(-22) ~ 2.7e-10
const double BMDS_MISSING = -9999.0;
const double BMDS_QNORM = 1.959964;  // bound for 95% confidence interval
extern std::string BMDS_VERSION;

enum nested_model { nlogistic = 1, nctr = 2 };

// BMDS helper structures
#ifdef _WIN64
#  pragma pack(8)
#elif _WIN32
#  pragma pack(4)
#endif

// forward declarations
double calcSlopeFactor(double bmr, double bmdl);

struct test_struct {
  double BMD;
  int n;
  bool validResult;
  std::vector<double> doses;
};
// BMD_results:
//   Purpose - Contains various BMD values returned by BMDS.
//   It is used to facilitate returning results needed for BMDS software.
struct BMDS_results {
  double BMD;
  double BMDL;
  double BMDU;
  double AIC;
  double BIC_equiv;
  double chisq;
  std::vector<bool> bounded;
  std::vector<double> stdErr;
  std::vector<double> lowerConf;
  std::vector<double> upperConf;
  bool validResult;
  double slopeFactor;

  void setSlopeFactor(double bmr) { slopeFactor = calcSlopeFactor(bmr, BMDL); }
};

struct BMDSMA_results {
  double BMD_MA;
  double BMDL_MA;
  double BMDU_MA;
  std::vector<double> BMD;
  std::vector<double> BMDL;
  std::vector<double> BMDU;
  std::vector<double> ebLower;  // size is number of dose groups
  std::vector<double> ebUpper;  // size is number of dose groups
};

// all arrays are length 4
struct testsOfInterest {
  std::vector<double> llRatio;
  std::vector<double> DF;
  std::vector<double> pVal;
};

// all arrays are of length 5
// Likelihoods of Interest
// indexing of arrays are:
// 0 - A1 Model
// 1 - A2 Model
// 2 - A3 Model
// 3 - fitted Model
// 4 - R Model
struct continuous_AOD {
  std::vector<double> LL;
  std::vector<int> nParms;
  std::vector<double> AIC;
  double addConst;
  struct testsOfInterest TOI;
};

struct dicho_AOD {
  double fullLL;  // A1;  //fullLL - Full Model
  int nFull;      // N1;     //NFull Full Model
  double redLL;   // A2;  //redLL Reduced Model
  int nRed;       // N2;     //NRed  Reduced Model
  double fittedLL;
  int nFit;
  double devFit;
  double devRed;
  int dfFit;
  int dfRed;
  double pvFit;
  double pvRed;
};

// each array has number of dose groups in suff_stat data
struct continuous_GOF {
  std::vector<double> dose;
  std::vector<double> size;
  std::vector<double> estMean;
  std::vector<double> calcMean;
  std::vector<double> obsMean;
  std::vector<double> estSD;
  std::vector<double> calcSD;
  std::vector<double> obsSD;
  std::vector<double> res;
  int n;  // total # of obs/doses
  std::vector<double> ebLower;
  std::vector<double> ebUpper;
};

struct dichotomous_GOF {
  int n;  // total number of observations obs/n
  std::vector<double> expected;
  std::vector<double> residual;
  double test_statistic;
  double p_value;
  double df;
  std::vector<double> ebLower;
  std::vector<double> ebUpper;
};

struct nestedBootstrap {
  std::vector<double> pVal;  // size = numRuns + 1
  std::vector<double> perc50;
  std::vector<double> perc90;
  std::vector<double> perc95;
  std::vector<double> perc99;
};

struct nestedLitterData {
  std::vector<double> dose;  // size = numRows
  std::vector<double> LSC;
  std::vector<double> estProb;
  std::vector<double> litterSize;
  std::vector<double> expected;
  std::vector<int> observed;
  std::vector<double> SR;
  double chiSq;
};

struct nestedReducedData {
  std::vector<double> dose;        // size = numRows
  std::vector<double> propAffect;  // estimate of proportion affected
  std::vector<double> lowerConf;
  std::vector<double> upperConf;
};

struct nestedSRData {
  double minSR;
  double avgSR;
  double maxSR;
  double minAbsSR;
  double avgAbsSR;
  double maxAbsSR;
};

struct python_dichotomous_analysis {
  int model;              // Model Type as listed in dich_model
  int n;                  // total number of observations obs/n
  std::vector<double> Y;  // observed +
  std::vector<double> doses;
  std::vector<double> n_group;  // size of the group
  std::vector<double> prior;    // a column order matrix (parms x prior_cols)
  int BMD_type;                 // 1 = extra ; added otherwise
  double BMR;
  double alpha;      // alpha of the analysis
  int degree;        // degree of polynomial used only  multistage
  int samples;       // number of MCMC samples.
  int burnin;        // size of burin
  int parms;         // number of parameters in the model
  int prior_cols;    // colunns in the prior
  bool penalizeAIC;  // whether to penalize the AIC by counting parameters that hit a bound
};

struct python_dichotomous_model_result {
  int model;                     // dichotomous model specification
  int nparms;                    // number of parameters in the model
  std::vector<double> parms;     // Parameter Estimate
  std::vector<double> cov;       // Covariance Estimate
  double max;                    // Value of the Likelihood/Posterior at the maximum
  int dist_numE;                 // number of entries in rows for the bmd_dist
  double model_df;               // Used model degrees of freedom
  double total_df;               // Total degrees of freedom
  std::vector<double> bmd_dist;  // bmd distribution (dist_numE x 2) matrix
  double bmd;                    // the central estimate of the BMD
  double gof_p_value;            // P-value from Chi Square goodness of fit
  double gof_chi_sqr_statistic;  // Chi Square Statistic for goodness of fit
  struct dichotomous_GOF gof;
  struct BMDS_results bmdsRes;
  struct dicho_AOD aod;

  double getSRAtDose(double targetDose, std::vector<double> doses);
};

struct python_dichotomousMA_analysis {
  int nmodels;                              // number of models for the model average
  std::vector<std::vector<double>> priors;  // List of pointers to prior arrays
                                            // priors[i] is the prior array for the ith model ect
  std::vector<int> nparms;                  // parameters in each model
  std::vector<int> actual_parms;            // actual number of parameters in the model
  std::vector<int> prior_cols;      // columns in the prior if there are 'more' in the future
                                    // presently there are only 5
  std::vector<int> models;          // list of models this is defined by dich_model.
  std::vector<double> modelPriors;  // prior probability on the model
  struct python_dichotomous_analysis pyDA;
};

struct python_dichotomousMA_result {
  int nmodels;  // number of models for each
  std::vector<python_dichotomous_model_result>
      models;                      // Individual model fits for each model average
  int dist_numE;                   // number of entries in rows for the bmd_dist
  std::vector<double> post_probs;  // posterior probabilities
  std::vector<double> bmd_dist;    // bmd ma distribution (dist_numE x 2) matrix
  struct BMDSMA_results bmdsRes;
};

struct python_continuous_analysis {
  enum cont_model model;
  int n;
  bool suff_stat;               // true if the data are in sufficient statistics format
  std::vector<double> Y;        // observed data means or actual data
  std::vector<double> doses;    //
  std::vector<double> sd;       // SD of the group if suff_stat = true, null otherwise.
  std::vector<double> n_group;  // N for each group if suff_stat = true, null otherwise
  std::vector<double> prior;    // a column order matrix px5 where p is the number of parametersd
  int BMD_type;                 // type of BMD
  bool isIncreasing;            // if the BMD is defined increasing or decreasing
  double BMR;                   // Benchmark response related to the BMD type
  double tail_prob;             // tail probability
  int disttype;                 // Distribution type defined in the enum distribution
  double alpha;                 // specified alpha
  int samples;                  // number of MCMC samples.
  int degree;                   // if polynomial it is the degree
  int burnin;                   // burn in
  int parms;                    // number of parameters
  int prior_cols;
  int transform_dose;  // Use the arc-sin-hyperbolic inverse to transform dose.
  bool restricted;
  bool detectAdvDir;
  bool penalizeAIC;  // whether to penalize the AIC by counting parameters that hit a bound
};

struct python_continuous_model_result {
  int model;                     // continuous model specification
  int dist;                      // distribution_type
  int nparms;                    // number of parameters in the model
  std::vector<double> parms;     // Parameter Estimate
  std::vector<double> cov;       // Covariance Estimate
  double max;                    // Value of the Likelihood/Posterior at the maximum
  int dist_numE;                 // number of entries in rows for the bmd_dist
  double model_df;               // Used model degrees of freedom
  double total_df;               // Total degrees of freedom
  double bmd;                    // The bmd at the maximum
  std::vector<double> bmd_dist;  // bmd distribution (dist_numE x 2) matrix
  struct continuous_GOF gof;
  struct BMDS_results bmdsRes;
  struct continuous_AOD aod;
};

struct python_multitumor_analysis {
  //  int model; // Model Type as listed in dich_model
  int ndatasets;

  std::vector<std::vector<python_dichotomous_analysis>> models;  //(size ndatasets * nmodels[i])

  std::vector<int> n;        // total number of observations per dataset (size ndatasets)
  std::vector<int> nmodels;  // # of models per dataset (size ndatasets)
  int BMD_type;              // 1 = extra ; added otherwise
  double BMR;
  double alpha;    // alpha of the analysis
  int prior_cols;  // colunns in the prior
  std::vector<int>
      degree;  // degree of selected polynomial used for each ind multistage (size ndatasets)
};

struct python_multitumor_result {
  int ndatasets;     // number of models for each
  bool validResult;  // true if BMD and slope factor were both calculated for multitumor model
  std::vector<int> nmodels;  // # of models per dataset (size ndatasets)
  std::vector<std::vector<python_dichotomous_model_result>>
      models;  // Individual model fits for each dataset nmodels[i]*ndatasets
  std::vector<int> selectedModelIndex;
  // int dist_numE; // number of entries in rows for the bmd_dist
  // std::vector<double> bmd_dist; // bmd ma distribution (dist_numE x 2) matrix
  double BMD;
  double BMDL;
  double BMDU;
  double slopeFactor;
  double combined_LL;        // combined log-likelihood
  double combined_LL_const;  // combined log-likelihood constant

  void setSlopeFactor(double bmr) { slopeFactor = calcSlopeFactor(bmr, BMDL); }
};

struct python_nested_analysis {
  enum nested_model model;  // model type in nest_model enum
  std::vector<double> doses;
  std::vector<double> litterSize;
  std::vector<double> incidence;
  std::vector<double> lsc;    // litter specific covariate
  std::vector<double> prior;  // a column order matrix (parms x prior_cols)
  int LSC_type;               // 1 = Overall Mean; 2 = control group mean; 0 = do not use LSC
  int ILC_type;               // 1 = estimate intralitter; assume 0 otherwise
  int BMD_type;               // 1 = extra;  added otherwise
  bool estBackground;         // if false, sets background to zero
  int parms;                  // number of parameters in model
  int prior_cols;
  double BMR;
  double alpha;
  int numBootRuns;   // number of bootstrap run
  int iterations;    // number of iterations per run
  long seed;         // -9999 = automatic;  seed value otherwise
  bool penalizeAIC;  // whether to penalize the AIC by counting parameters that hit a bound
};

struct python_nested_result {
  bool validResult;
  enum nested_model model;
  int nparms;
  std::vector<double> parms;
  std::vector<double> cov;
  int dist_numE;    // number of entries in rows for the bmd dist
  double model_df;  // Used model degrees of freedom
  double bmd;
  double fixedLSC;
  double LL;
  double combPVal;
  struct BMDS_results bmdsRes;
  // NestedLitterDataRow
  struct nestedLitterData litter;
  // NestedBootstrapRow
  struct nestedBootstrap boot;
  // Nested Reduced DataRow
  struct nestedReducedData reduced;
  // Nested Scaled Residual Data
  struct nestedSRData srData;
};

struct msComboEq {
  double bmr;
  int nT;
  std::vector<int> degree;
};

struct msComboInEq {
  int nT;
  double target;
  std::vector<int> nObs;
  std::vector<int> degree;
  std::vector<std::vector<double>> doses;
  std::vector<std::vector<double>> Y;
  std::vector<std::vector<double>> n_group;
};

// for nested models
struct nested_AOD {
  double LL;
  double dev;
  double df;
  double pv;
};

struct nestedObjData {
  std::vector<double> Ls;   // Litter size
  std::vector<double> Xi;   // Dose
  std::vector<int> Xg;      // Dose group
  std::vector<double> Yp;   // positive response
  std::vector<double> Yn;   // negative response
  std::vector<double> Lsc;  // Litter specific covariate
  std::vector<double> prior;
  std::vector<double> GXi;  // doses at each dose group
  std::vector<bool> Spec;
  int ngrp;
  double smax;
  double smin;
  double smean;
  int LSC_type;
  int ILC_type;
  double isBMDL;
  // only used for BMDL
  double ck;
  double LR;
  double xlk;     // tmp likelihood (for BMDL calc)
  double BMD_lk;  // likelihood for BMD
  double tD;      // initially holds BMD dose
  double sijfixed;
  int riskType;
  double BMR;
  double tol;     // tolerance for optimization
  int optimizer;  // 1=LD_SLSQP, 2=???, 3=LD_LBFGS
  enum nested_model model;
};

#ifdef _WIN32
#  pragma pack()
#endif

// c entry
#ifdef __cplusplus
extern "C" {
#endif

void cleanDouble(double *val);

void rescale_dichoParms(int model, double *parms);
void rescale_contParms(struct continuous_analysis *CA, double *parms);

void calcParmCIs_dicho(struct dichotomous_model_result *res, struct BMDS_results *bmdsRes);
void calcParmCIs_cont(struct continuous_model_result *res, struct BMDS_results *bmdsRes);

void bmdsConvertSStat(
    struct continuous_analysis *ca, struct continuous_analysis *newCA, bool clean
);

void calcTestsOfInterest(struct continuous_AOD *aod);

void determineAdvDir(struct continuous_analysis *anal);

void calc_contAOD(
    struct continuous_analysis *CA, struct continuous_analysis *GOFanal,
    struct continuous_model_result *res, struct BMDS_results *bmdsRes, struct continuous_AOD *aod
);
void calc_dichoAOD(
    struct dichotomous_analysis *DA, struct dichotomous_model_result *res,
    struct BMDS_results *bmdsRes, struct dicho_AOD *bmdsAOD, struct dichotomous_aod *aod
);

void collect_dicho_bmd_values(
    struct dichotomous_analysis *anal, struct dichotomous_model_result *res,
    struct BMDS_results *BMDres
);
void collect_dichoMA_bmd_values(
    struct dichotomousMA_analysis *anal, struct dichotomousMA_result *res,
    struct BMDSMA_results *BMDres, double alpha
);
void collect_cont_bmd_values(
    struct continuous_analysis *anal, struct continuous_model_result *res,
    struct BMDS_results *BMDres
);

void calcDichoAIC(
    struct dichotomous_analysis *anal, struct dichotomous_model_result *res,
    struct BMDS_results *BMDSres, bool penalizeAIC
);

void calcContAIC(
    struct continuous_analysis *anal, struct continuous_model_result *res,
    struct BMDS_results *BMDSres, bool penalizeAIC
);

double calcNestedAIC(
    double fitted_LL, double fitted_df, double red_df, int numBounded, bool penalizeAIC
);

void clean_dicho_results(
    struct dichotomous_model_result *res, struct dichotomous_GOF *gof, struct BMDS_results *bmdsRes,
    struct dicho_AOD *aod
);
void clean_cont_results(
    struct continuous_model_result *res, struct BMDS_results *bmdsRes, struct continuous_AOD *aod,
    struct continuous_GOF *gof
);
void clean_dicho_MA_results(struct dichotomousMA_result *res, struct BMDSMA_results *bmdsRes);

void clean_multitumor_results(struct python_multitumor_result *res);

void clean_nested_results(struct python_nested_result *res);

void convertFromPythonDichoAnalysis(
    struct dichotomous_analysis *anal, struct python_dichotomous_analysis *pyAnal
);
void convertToPythonDichoRes(
    struct dichotomous_model_result *res, struct python_dichotomous_model_result *pyRes
);
void convertFromPythonDichoRes(
    struct dichotomous_model_result *res, struct python_dichotomous_model_result *ret
);

void selectMultitumorModel(
    struct python_multitumor_analysis *pyAnal, struct python_multitumor_result *pyRes
);
int selectBestMultitumorModel(
    std::vector<python_dichotomous_analysis> &analModels,
    std::vector<python_dichotomous_model_result> &resModels
);

void BMDS_ENTRY_API __stdcall runMultitumorModel(
    struct python_multitumor_analysis *pyAnal, struct python_multitumor_result *pyRes
);

double DLgamma(double x);

double LogLik_Constant(std::vector<double> Y, std::vector<double> n_group);

double zeroin(
    double ax, double bx, double tol, double (*f)(int, double[], double, double), int nparm,
    double Parms[], double ck
);
double zeroin_nested(
    double ax, double bx, double tol,
    double (*f)(std::vector<double> &, double, double, struct nestedObjData *),
    std::vector<double> &Parms, double ck, struct nestedObjData *objData
);

double BMD_func(int n, double p[], double x, double ck);

double getclmt(
    python_multitumor_analysis *pyAnal, python_multitumor_result *pyRes, double Dose, double target,
    double maxDose, std::vector<double> xParms, std::vector<double> fixedParm, bool isBMDL
);

void multitumorCLs(
    struct python_multitumor_analysis *pyAnal, struct python_multitumor_result *pyRes, double Dose,
    double D, double LR, double gtol, int *is_zero
);

void Multistage_ComboBMD(
    struct python_multitumor_analysis *pyAnal, struct python_multitumor_result *pyRes
);

double objfunc_bmdl(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data);
double objfunc_bmdu(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data);
double objfunc2(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data);
double myEqualityConstraint(const std::vector<double> &x, std::vector<double> &grad, void *data);
double myInequalityConstraint1(const std::vector<double> &x, std::vector<double> &grad, void *data);
double myInequalityConstraint2(const std::vector<double> &x, std::vector<double> &grad, void *data);
double myInequalityConstraint3(const std::vector<double> &x, std::vector<double> &grad, void *data);

double slog(double X);
double dslog(double P);

double round_to(double value, double precision = 1.0);

double ComboMaxLike2(
    int flag, double dose, double *crisk, std::vector<std::vector<double>> p,
    python_multitumor_analysis *pyAnal, python_multitumor_result *pyres
);

void validatePositiveInput(double &val);

void probability_inrange(double *ex);

double Nlogist_lk(std::vector<double> p, struct nestedObjData *objData);
double NCTR_lk(std::vector<double> p, struct nestedObjData *objData);

double Nlogist_g(std::vector<double> p, std::vector<double> &g, struct nestedObjData *objData);
double NCTR_g(std::vector<double> p, std::vector<double> &g, struct nestedObjData *objData);

void Nlogist_probs(
    std::vector<double> &probs, const std::vector<double> &p, bool compgrad,
    std::vector<std::vector<double>> &gradij, struct nestedObjData *objData
);
void NCTR_probs(
    std::vector<double> &probs, const std::vector<double> &p, bool compgrad,
    std::vector<std::vector<double>> &gradij, struct nestedObjData *objData
);

double opt_nlogistic(std::vector<double> &p, struct nestedObjData *data);
double opt_nctr(std::vector<double> &p, struct nestedObjData *data);

double objfunc_nlogistic_ll(const std::vector<double> &p, std::vector<double> &grad, void *data);
double objfunc_nctr_ll(const std::vector<double> &p, std::vector<double> &grad, void *data);

double nestedInequalityConstraint(
    const std::vector<double> &x, std::vector<double> &grad, void *data
);

void Nlogist_BMD(
    struct python_nested_analysis *pyAnal, struct python_nested_result *pyRes, double smin,
    double smax, double sijfixed, double xmax, struct nestedObjData *objData
);
void Nctr_BMD(
    struct python_nested_analysis *pyAnal, struct python_nested_result *pyRes, double smin,
    double smax, double sijfixed, double xmax, struct nestedObjData *objData
);

void Nlogist_vcv(
    std::vector<double> &p, std::vector<bool> &bounded, struct nestedObjData *objData,
    std::vector<std::vector<double>> &vcv
);
void NCTR_vcv(
    std::vector<double> &p, std::vector<bool> &bounded, struct nestedObjData *objData,
    std::vector<std::vector<double>> &vcv
);

void Nlogist_grad(std::vector<double> &p, struct nestedObjData *objData, std::vector<double> &grad);
void NCTR_grad(std::vector<double> &p, struct nestedObjData *objData, std::vector<double> &grad);

double BMDL_func_nlog(std::vector<double> &p, double D, double gtol, struct nestedObjData *objData);
double BMDL_func_nctr(std::vector<double> &p, double D, double gtol, struct nestedObjData *objData);

double QCHISQ(double p, int m);
double CHISQ(double x, int m);

void outputObjData(struct nestedObjData *objData);

void Nested_Bootstrap(
    struct nestedObjData *objData, struct python_nested_analysis *pyAnal,
    struct python_nested_result *pyRes, long seed, int iterations, int BSLoops
);

void Nested_SRoI(
    struct nestedObjData *objData, struct nestedSRData *srData, std::vector<double> &SR,
    const std::vector<int> &grpSize, double bmd
);

void Nested_GOF(
    const std::vector<double> &parms, struct nestedObjData *objData,
    struct nestedLitterData *litterData, const std::vector<int> &grpSize
);

void Nested_reduced(double alpha, struct nestedObjData *objData, struct nestedReducedData *redData);

void raoscott(struct nestedObjData *objData, std::vector<double> &num, std::vector<double> &den);

void Quantal_CI(
    std::vector<double> &Yp, std::vector<double> &Yn, double conf, struct nestedReducedData *redData
);

void Nlogist_Predict(
    const std::vector<double> &parms, struct nestedObjData *objData, std::vector<double> &P
);
void NCTR_Predict(
    const std::vector<double> &parms, struct nestedObjData *objData, std::vector<double> &P
);

double calcNlogisticCLs(
    double xa, double xb, std::vector<double> &pint, struct nestedObjData *objData, bool isLower
);
double calcNCTRCLs(
    double xa, double xb, std::vector<double> &pint, struct nestedObjData *objData, bool isLower
);

void SortNestedData(
    const std::vector<int> &GrpSize, std::vector<double> &Xi, std::vector<double> &Ls,
    std::vector<double> &Yp, std::vector<double> &Yn, std::vector<double> &Lsc, bool sortByLsc
);

void BMDS_ENTRY_API __stdcall runBMDSDichoAnalysis(
    struct dichotomous_analysis *anal, struct dichotomous_model_result *res,
    struct dichotomous_GOF *gof, struct BMDS_results *bmdsRes, struct dicho_AOD *aod,
    bool *penalizeAIC
);

void BMDS_ENTRY_API __stdcall runBMDSContAnalysis(
    struct continuous_analysis *anal, struct continuous_model_result *res,
    struct BMDS_results *bmdsRes, struct continuous_AOD *aod, struct continuous_GOF *gof,
    bool *detectAdvDir, bool *restricted, bool *penalizeAIC
);

void BMDS_ENTRY_API __stdcall runBMDSDichoMA(
    struct dichotomousMA_analysis *MA, struct dichotomous_analysis *DA,
    struct dichotomousMA_result *res, struct BMDSMA_results *bmdsRes
);

string BMDS_ENTRY_API __stdcall version();

void BMDS_ENTRY_API __stdcall pythonBMDSDicho(
    struct python_dichotomous_analysis *pyAnal, struct python_dichotomous_model_result *pyRes
);

void BMDS_ENTRY_API __stdcall pythonBMDSDichoMA(
    struct python_dichotomousMA_analysis *pyMA, struct python_dichotomousMA_result *pyRes
);

void BMDS_ENTRY_API __stdcall pythonBMDSCont(
    struct python_continuous_analysis *pyAnal, struct python_continuous_model_result *pyRes
);

void BMDS_ENTRY_API __stdcall pythonBMDSMultitumor(
    struct python_multitumor_analysis *pyAnal, struct python_multitumor_result *pyRes
);

void BMDS_ENTRY_API __stdcall pythonBMDSNested(
    struct python_nested_analysis *pyAnal, struct python_nested_result *pyRes
);

#ifdef __cplusplus
}
#endif

// overloaded print statements
std::string BMDS_ENTRY_API __stdcall printBmdsStruct(
    struct BMDS_results *bmdsRes, bool print = true
);

std::string BMDS_ENTRY_API __stdcall printBmdsStruct(
    struct BMDSMA_results *bmdsRes, bool print = true
);

std::string BMDS_ENTRY_API __stdcall printBmdsStruct(
    struct nestedLitterData *litter, bool print = true
);

std::string BMDS_ENTRY_API __stdcall printBmdsStruct(
    struct nestedBootstrap *boot, bool print = true
);

std::string BMDS_ENTRY_API __stdcall printBmdsStruct(
    struct nestedReducedData *red, bool print = true
);

std::string BMDS_ENTRY_API __stdcall printBmdsStruct(struct nestedSRData *sr, bool print = true);

std::string BMDS_ENTRY_API __stdcall printBmdsStruct(
    struct python_continuous_analysis *pyAnal, bool print = true
);

std::string BMDS_ENTRY_API __stdcall printBmdsStruct(
    struct python_continuous_model_result *pyRes, bool print = true
);

std::string BMDS_ENTRY_API __stdcall printBmdsStruct(
    struct python_dichotomousMA_analysis *pyMA, bool print = true
);

std::string BMDS_ENTRY_API __stdcall printBmdsStruct(
    struct python_dichotomousMA_result *pyRes, bool print = true
);

std::string BMDS_ENTRY_API __stdcall printBmdsStruct(
    struct python_multitumor_analysis *pyAnal, bool print = true
);

std::string BMDS_ENTRY_API __stdcall printBmdsStruct(
    struct python_multitumor_result *pyRes, bool print = true
);

std::string BMDS_ENTRY_API __stdcall printBmdsStruct(
    struct python_nested_analysis *pyAnal, bool print = true
);

std::string BMDS_ENTRY_API __stdcall printBmdsStruct(
    struct python_nested_result *pyRes, bool print = true
);

std::string BMDS_ENTRY_API __stdcall printBmdsStruct(struct continuous_GOF *gof, bool print = true);

std::string BMDS_ENTRY_API __stdcall printBmdsStruct(struct continuous_AOD *AOD, bool print = true);

std::string BMDS_ENTRY_API __stdcall printBmdsStruct(
    struct testsOfInterest *TOI, bool print = true
);

std::string BMDS_ENTRY_API __stdcall printBmdsStruct(
    struct python_dichotomous_analysis *pyAnal, bool print = true
);

std::string BMDS_ENTRY_API __stdcall printBmdsStruct(
    struct python_dichotomous_model_result *pyRes, bool print = true
);

std::string BMDS_ENTRY_API __stdcall printBmdsStruct(
    struct dichotomous_GOF *gof, bool print = true
);

std::string BMDS_ENTRY_API __stdcall printBmdsStruct(struct dicho_AOD *AOD, bool print = true);
