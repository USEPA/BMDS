#ifdef __cplusplus
  #include <Eigen/Dense>
#endif
//#include "bmdStruct.h"
#include "dichotomous_entry_code.h"

const double BMDS_EPS = 1.0e-6;

// BMDS helper structures

// BMD_results:
//   Purpose - Contains various BMD values returned by BMDS.
//   It is used to facilitate returning results needed for BMDS software. 
struct BMDS_results{
  double BMD;   
  double BMDL; 
  double BMDU; 
  double AIC;   
  bool *bounded;
};


//enum TestOfInterest {A1 = 1, A2 = 2, A3 = 3, fitted = 4, R = 5};

struct testsOfInterest {
  double *llRatio;
  double *DF;
  double *pVal;
};
struct continuous_AOD{
  double *LL;
  double *nParms;
  double *AIC;
  double addConst;
  struct testsOfInterest *TOI;
};


struct testsOfInterest1{
  double *logLikeRatio;
  double *tmp;
//  double *DF;
//  double *pValue;
};

//
struct continuous_AOD1{
  double *logLikelihood;
  double *nparms;
  double *AIC;
  double addConst;
  //struct testsOfInterest TOI;
};


//c entry
#ifdef __cplusplus
extern "C" {
#endif

void bmdsConvertSStat(struct continuous_analysis *ca, struct continuous_analysis *newCA);

void calcTestsOfInterest(struct continuous_AOD *aod);

void determineAdvDir(struct continuous_analysis *anal);

void calc_AOD(struct continuous_analysis *CA, struct continuous_model_result *res, struct BMDS_results *bmdsRes, struct continuous_AOD *aod);

//void collect_dicho_bmd_values(double *bmd_dist, struct BMD_results *BMDres);
void collect_dicho_bmd_values(struct dichotomous_analysis *anal, struct dichotomous_model_result *res, struct BMDS_results *BMDres);

void collect_cont_bmd_values(struct continuous_analysis *anal, struct continuous_model_result *res, struct BMDS_results *BMDres);


void runBMDSDichoAnalysis(struct dichotomous_analysis *anal, struct dichotomous_model_result *res, struct dichotomous_PGOF_result *gofRes, struct BMDS_results *bmdsRes);


void runBMDSContAnalysis(struct continuous_analysis *anal, struct continuous_model_result *res, struct BMDS_results *bmdsRes, struct continuous_AOD *aod, bool detectAdvDir);
#ifdef __cplusplus
}
#endif
