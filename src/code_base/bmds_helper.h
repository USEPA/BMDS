#ifdef __cplusplus
  #include <Eigen/Dense>
#endif
//#include "bmdStruct.h"
#include "dichotomous_entry_code.h"

//define __stdcall to nothing if not on Windows platform
#ifndef WIN32
#define __stdcall
#endif

//define import/export attributes if on Windows platform
#ifndef R_COMPILATION
#  ifndef _WIN32
#    define BMDS_ENTRY_API
#  else
#    ifdef RBMDS_EXPORTS
#      define BMDS_ENTRY_API __declspec(dllexport)
#    else
#      define BMDS_ENTRY_API __declspec(dllimport)
#    endif // BMDS_MODELS_EXPORTS
#  endif
#else
#  define BMDS_ENTRY_API
#endif

const double BMDS_EPS = 1.0e-6;
const double BMDS_MISSING = -9999.0;
const double BMDS_QNORM = 1.959964;  //bound for 95% confidence interval
extern std::string BMDS_VERSION; 

// BMDS helper structures
#ifdef _WIN64
#pragma pack(8)
#elif _WIN32
#pragma pack(4)
#endif

struct test_struct{
  double BMD;
  int n;
  bool validResult;
  std::vector<double> doses;
};
// BMD_results:
//   Purpose - Contains various BMD values returned by BMDS.
//   It is used to facilitate returning results needed for BMDS software. 
struct BMDS_results{
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
};

struct BMDSMA_results{
  double BMD_MA;
  double BMDL_MA;
  double BMDU_MA;
  std::vector<double> BMD;
  std::vector<double> BMDL;
  std::vector<double> BMDU;
  std::vector<double> ebLower;  //size is number of dose groups
  std::vector<double> ebUpper;  //size is number of dose groups
};


//all arrays are length 4
struct testsOfInterest {
  std::vector<double> llRatio;
  std::vector<double> DF;
  std::vector<double> pVal;
};


//all arrays are of length 5
//Likelihoods of Interest
//indexing of arrays are:
//0 - A1 Model
//1 - A2 Model
//2 - A3 Model
//3 - fitted Model
//4 - R Model
struct continuous_AOD{
  std::vector<double> LL;
  std::vector<int> nParms;
  std::vector<double> AIC;
  double addConst;
  struct testsOfInterest *TOI;
};

struct dicho_AOD{
  double fullLL;  //A1;  //fullLL - Full Model
  int nFull;      //N1;     //NFull Full Model
  double redLL;   //A2;  //redLL Reduced Model
  int nRed;       //N2;     //NRed  Reduced Model
  double fittedLL;
  int nFit;
  double devFit;
  double devRed;
  int dfFit;
  int dfRed;
  double pvFit;
  double pvRed;
};


//each array has number of dose groups in suff_stat data
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
  int n; //total # of obs/doses  
  std::vector<double> ebLower;
  std::vector<double> ebUpper;
};

struct dichotomous_GOF {
  int     n;        // total number of observations obs/n 
  std::vector<double> expected; 
  std::vector<double> residual;
  double  test_statistic; 
  double  p_value; 
  double  df;  
  std::vector<double> ebLower;
  std::vector<double> ebUpper;
};

struct python_dichotomous_analysis{
  int model; // Model Type as listed in dich_model
  int n;     // total number of observations obs/n
  std::vector<double> Y;  // observed +
  std::vector<double> doses;  
  std::vector<double> n_group; //size of the group
  std::vector<double> prior;  //a column order matrix (parms x prior_cols)
  int BMD_type; // 1 = extra ; added otherwise
  double BMR;
  double alpha; // alpha of the analysis
  int degree;  // degree of polynomial used only  multistage
  int samples; // number of MCMC samples.
  int burnin;  // size of burin
  int parms;   // number of parameters in the model
  int prior_cols; // colunns in the prior
};

struct python_dichotomous_model_result{
  int      model;               // dichotomous model specification
  int      nparms;              //number of parameters in the model
  std::vector<double> parms;    // Parameter Estimate
  std::vector<double> cov;      // Covariance Estimate
  double   max;                 // Value of the Likelihood/Posterior at the maximum
  int      dist_numE;           // number of entries in rows for the bmd_dist
  double      model_df;         // Used model degrees of freedom
  double      total_df;         // Total degrees of freedom
  std::vector<double> bmd_dist; // bmd distribution (dist_numE x 2) matrix
  double  bmd;                  // the central estimate of the BMD
  double gof_p_value;           // P-value from Chi Square goodness of fit
  double gof_chi_sqr_statistic; // Chi Square Statistic for goodness of fit
  struct dichotomous_GOF gof;
  struct BMDS_results bmdsRes;
  struct dicho_AOD aod;
};

struct python_dichotomousMA_analysis{
  int    nmodels;      		//number of models for the model average
  std::vector<std::vector<double>> priors;     		// List of pointers to prior arrays
                       		// priors[i] is the prior array for the ith model ect
  std::vector<int> nparms;     //parameters in each model
  std::vector<int> actual_parms;//actual number of parameters in the model
  std::vector<int> prior_cols;  // columns in the prior if there are 'more' in the future
                       		// presently there are only 5
  std::vector<int> models;      // list of models this is defined by dich_model.
  std::vector<double> modelPriors; // prior probability on the model
  struct python_dichotomous_analysis pyDA;
};

struct python_dichotomousMA_result{
  int                       nmodels; //number of models for each
  std::vector<python_dichotomous_model_result> models;  //Individual model fits for each model average
  int dist_numE; // number of entries in rows for the bmd_dist
  std::vector<double> post_probs; // posterior probabilities
  std::vector<double> bmd_dist; // bmd ma distribution (dist_numE x 2) matrix
  struct BMDSMA_results bmdsRes;
};

struct python_continuous_analysis{
  enum cont_model model;
  int n;
  bool suff_stat; //true if the data are in sufficient statistics format
  std::vector<double> Y; // observed data means or actual data
  std::vector<double> doses; //
  std::vector<double> sd; // SD of the group if suff_stat = true, null otherwise.
  std::vector<double> n_group; // N for each group if suff_stat = true, null otherwise
  std::vector<double> prior; // a column order matrix px5 where p is the number of parametersd
  int BMD_type; // type of BMD
  bool isIncreasing; // if the BMD is defined increasing or decreasing
  double BMR; // Benchmark response related to the BMD type
  double tail_prob; // tail probability
  int    disttype;  // Distribution type defined in the enum distribution
  double alpha;     // specified alpha
  int samples; // number of MCMC samples.
  int degree; // if polynomial it is the degree
  int burnin;  // burn in
  int parms; // number of parameters
  int prior_cols;
  int transform_dose; // Use the arc-sin-hyperbolic inverse to transform dose.
  bool restricted;
  bool detectAdvDir;
};

struct python_continuous_model_result{
  int      model;           // continuous model specification
  int      dist;            // distribution_type
  int      nparms;                    //number of parameters in the model
  std::vector<double>  parms;           // Parameter Estimate
  std::vector<double>  cov;             // Covariance Estimate
  double   max;             // Value of the Likelihood/Posterior at the maximum
  int      dist_numE;       // number of entries in rows for the bmd_dist
  double    model_df;        // Used model degrees of freedom
  double    total_df;        // Total degrees of freedom
  double    bmd;             // The bmd at the maximum
  std::vector<double> bmd_dist;        // bmd distribution (dist_numE x 2) matrix
  struct continuous_GOF gof;
  struct BMDS_results bmdsRes;
  struct continuous_AOD aod;
};


struct python_multitumor_analysis{
//  int model; // Model Type as listed in dich_model
  int ndatasets;

  std::vector<std::vector<python_dichotomous_analysis>> models; //(size ndatasets * nmodels[i])

  std::vector<int> n;     // total number of observations per dataset (size ndatasets)
  std::vector<int> nmodels;  //# of models per dataset (size ndatasets)
  //std::vector<std::vector<double>> Y;  //observed + (size ndatasets*n[i])
  //std::vector<std::vector<double>> doses;  //size ndatasets*n[i]
  //std::vector<std::vector<double>> n_group; //size of group (size ndatasets*n[i])
  int BMD_type; // 1 = extra ; added otherwise
  double BMR;
  double alpha; // alpha of the analysis
  int prior_cols; // colunns in the prior
  std::vector<int> degree;  // degree of selected polynomial used for each ind multistage (size ndatasets)
  //std::vector<int> parms;   // number of parameters in the model
  std::vector<std::vector<double>> prior;  //a column order matrix (parms x prior_cols)
};

struct python_multitumor_result{
  int ndatasets; //number of models for each
  std::vector<int> nmodels; //# of models per dataset (size ndatasets)
  std::vector<std::vector<python_dichotomous_model_result>> models;  //Individual model fits for each dataset nmodels[i]*ndatasets
  std::vector<std::vector<dichotomous_GOF>> gofs;
  std::vector<std::vector<BMDS_results>> bmdsRess;
  std::vector<std::vector<dicho_AOD>> aods;  
  int dist_numE; // number of entries in rows for the bmd_dist
  std::vector<double> bmd_dist; // bmd ma distribution (dist_numE x 2) matrix
};

struct BMDSmultitumor_results{
  double BMD_MA;
  double BMDL_MA;
  double BMDU_MA;
  int nFitted;  //# of fitted datasets (may be less than number of submitted datasets)
  std::vector<double> BMD;  //size nFitted
  std::vector<double> BMDL; //size nFitted
  std::vector<double> BMDU; //size nFitted
  double combined_LL;  //combined log-likelihood 
  double combined_LL_const; //combined log-likelihood constant
};



#ifdef _WIN32
#pragma pack()
#endif

//c entry
#ifdef __cplusplus
extern "C" {
#endif

void cleanDouble(double *val);

void rescale_dichoParms(int model, double *parms);

void rescale_contParms(struct continuous_analysis *CA, double *parms);

void calcParmCIs_dicho(struct dichotomous_model_result *res, struct BMDS_results *bmdsRes);

void calcParmCIs_cont(struct continuous_model_result *res, struct BMDS_results *bmdsRes);

void bmdsConvertSStat(struct continuous_analysis *ca, struct continuous_analysis *newCA, bool clean);

void calcTestsOfInterest(struct continuous_AOD *aod);

void determineAdvDir(struct continuous_analysis *anal);

void calc_contAOD(struct continuous_analysis *CA, struct continuous_analysis *GOFanal,  struct continuous_model_result *res, struct BMDS_results *bmdsRes, struct continuous_AOD *aod);

void calc_dichoAOD(struct dichotomous_analysis *DA, struct dichotomous_model_result *res, struct BMDS_results *bmdsRes, struct dicho_AOD *bmdsAOD, struct dichotomous_aod *aod);

void collect_dicho_bmd_values(struct dichotomous_analysis *anal, struct dichotomous_model_result *res, struct BMDS_results *BMDres, double estParmCount);

void collect_dichoMA_bmd_values(struct dichotomousMA_analysis *anal, struct dichotomousMA_result *res, struct BMDSMA_results *BMDres);

void collect_cont_bmd_values(struct continuous_analysis *anal, struct continuous_model_result *res, struct BMDS_results *BMDres);

void clean_dicho_results(struct dichotomous_model_result *res, struct dichotomous_GOF *gof, struct BMDS_results *bmdsRes, struct dicho_AOD *aod);

void clean_cont_results(struct continuous_model_result *res, struct BMDS_results *bmdsRes, struct continuous_AOD *aod, struct continuous_GOF *gof);

void clean_dicho_MA_results(struct dichotomousMA_result *res, struct BMDSMA_results *bmdsRes);


void convertFromPythonDichoAnalysis(struct dichotomous_analysis *anal, struct python_dichotomous_analysis *pyAnal);

void convertToPythonDichoRes(struct dichotomous_model_result *res, struct python_dichotomous_model_result *pyRes);

void convertFromPythonDichoRes(struct dichotomous_model_result *res, struct python_dichotomous_model_result *ret);

void BMDS_ENTRY_API __stdcall runBMDSDichoAnalysis(struct dichotomous_analysis *anal, struct dichotomous_model_result *res, struct dichotomous_GOF *gof, struct BMDS_results *bmdsRes, struct dicho_AOD *aod);

void BMDS_ENTRY_API __stdcall runBMDSContAnalysis(struct continuous_analysis *anal, struct continuous_model_result *res, struct BMDS_results *bmdsRes, struct continuous_AOD *aod, struct continuous_GOF *gof, bool *detectAdvDir, bool *restricted);

void BMDS_ENTRY_API __stdcall runBMDSDichoMA(struct dichotomousMA_analysis *MA, struct dichotomous_analysis *DA,  struct dichotomousMA_result *res, struct BMDSMA_results *bmdsRes);

string BMDS_ENTRY_API __stdcall version();

int BMDS_ENTRY_API __stdcall add2(int i, int j);

void BMDS_ENTRY_API __stdcall testFun(struct test_struct *t);

void BMDS_ENTRY_API __stdcall pythonBMDSDicho(struct python_dichotomous_analysis *pyAnal, struct python_dichotomous_model_result *pyRes);

void BMDS_ENTRY_API __stdcall pythonBMDSDichoMA(struct python_dichotomousMA_analysis *pyMA, struct python_dichotomousMA_result *pyRes);

void BMDS_ENTRY_API __stdcall pythonBMDSCont(struct python_continuous_analysis *pyAnal, struct python_continuous_model_result *pyRes);

void BMDS_ENTRY_API __stdcall pythonBMDSMultitumor(struct python_multitumor_analysis *pyAnal, struct python_multitumor_result *pyRes, struct BMDSmultitumor_results *bmdsRes);

#ifdef __cplusplus
}
#endif
