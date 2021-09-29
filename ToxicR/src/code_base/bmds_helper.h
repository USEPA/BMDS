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
const char BMDS_VERSION[32] = "2021.09"; 

// BMDS helper structures
#ifdef _WIN64
#pragma pack(8)
#elif _WIN32
#pragma pack(4)
#endif
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
  bool *bounded;
  double *stdErr;
  double *lowerConf;
  double *upperConf;
  bool validResult;
};

struct BMDSMA_results{
  double BMD_MA;
  double BMDL_MA;
  double BMDU_MA;
  double *BMD;
  double *BMDL;
  double *BMDU;
  double *ebLower;  //size is number of dose groups
  double *ebUpper;  //size is number of dose groups
};


//all arrays are length 4
struct testsOfInterest {
  double *llRatio;
  double *DF;
  double *pVal;
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
  double *LL;
  int *nParms;
  double *AIC;
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
  double *dose;
  double *size;
  double *estMean;
  double *calcMean;
  double *obsMean;
  double *estSD;
  double *calcSD;
  double *obsSD;
  double *res;
  int n; //total # of obs/doses  
  double *ebLower;
  double *ebUpper;
};

struct dichotomous_GOF {
  int     n;        // total number of observations obs/n 
  double *expected; // 
  double *residual; //size of the group
  double  test_statistic; 
  double  p_value; 
  double  df;  
  double *ebLower;
  double *ebUpper;
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

void bmdsConvertSStat(struct continuous_analysis *ca, struct continuous_analysis *newCA);

void calcTestsOfInterest(struct continuous_AOD *aod);

void determineAdvDir(struct continuous_analysis *anal);

void calc_contAOD(struct continuous_analysis *CA, struct continuous_model_result *res, struct BMDS_results *bmdsRes, struct continuous_AOD *aod);

void calc_dichoAOD(struct dichotomous_analysis *DA, struct dichotomous_model_result *res, struct BMDS_results *bmdsRes, struct dicho_AOD *bmdsAOD);

//void collect_dicho_bmd_values(double *bmd_dist, struct BMD_results *BMDres);
void collect_dicho_bmd_values(struct dichotomous_analysis *anal, struct dichotomous_model_result *res, struct BMDS_results *BMDres, double estParmCount);

void collect_dichoMA_bmd_values(struct dichotomousMA_analysis *anal, struct dichotomousMA_result *res, struct BMDSMA_results *BMDres);

void collect_cont_bmd_values(struct continuous_analysis *anal, struct continuous_model_result *res, struct BMDS_results *BMDres);

void clean_dicho_results(struct dichotomous_model_result *res, struct dichotomous_GOF *gof, struct BMDS_results *bmdsRes, struct dicho_AOD *aod);

void clean_cont_results(struct continuous_model_result *res, struct BMDS_results *bmdsRes, struct continuous_AOD *aod, struct continuous_GOF *gof);

void clean_dicho_MA_results(struct dichotomousMA_result *res, struct BMDSMA_results *bmdsRes);

void BMDS_ENTRY_API __stdcall runBMDSDichoAnalysis(struct dichotomous_analysis *anal, struct dichotomous_model_result *res, struct dichotomous_GOF *gof, struct BMDS_results *bmdsRes, struct dicho_AOD *aod);

void BMDS_ENTRY_API __stdcall runBMDSContAnalysis(struct continuous_analysis *anal, struct continuous_model_result *res, struct BMDS_results *bmdsRes, struct continuous_AOD *aod, struct continuous_GOF *gof, bool *detectAdvDir, bool *restricted);

void BMDS_ENTRY_API __stdcall runBMDSDichoMA(struct dichotomousMA_analysis *MA, struct dichotomous_analysis *DA,  struct dichotomousMA_result *res, struct BMDSMA_results *bmdsRes);

void BMDS_ENTRY_API __stdcall version(char * versionStr);


#ifdef __cplusplus
}
#endif
