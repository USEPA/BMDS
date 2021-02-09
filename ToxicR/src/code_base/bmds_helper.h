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
  double chisq;
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
struct dicho_AOD{
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
  double pvFit;
  double pvRed;
};

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
};

//c entry
#ifdef __cplusplus
extern "C" {
#endif

void bmdsConvertSStat(struct continuous_analysis *ca, struct continuous_analysis *newCA);

void calcTestsOfInterest(struct continuous_AOD *aod);

void determineAdvDir(struct continuous_analysis *anal);

void calc_contAOD(struct continuous_analysis *CA, struct continuous_model_result *res, struct BMDS_results *bmdsRes, struct continuous_AOD *aod);

void calc_dichoAOD(struct dichotomous_analysis *DA, struct dichotomous_model_result *res, struct BMDS_results *bmdsRes, struct dicho_AOD *bmdsAOD);

//void collect_dicho_bmd_values(double *bmd_dist, struct BMD_results *BMDres);
void collect_dicho_bmd_values(struct dichotomous_analysis *anal, struct dichotomous_model_result *res, struct BMDS_results *BMDres);

void collect_cont_bmd_values(struct continuous_analysis *anal, struct continuous_model_result *res, struct BMDS_results *BMDres);


void runBMDSDichoAnalysis(struct dichotomous_analysis *anal, struct dichotomous_model_result *res, struct dichotomous_PGOF_result *gofRes, struct BMDS_results *bmdsRes, struct dicho_AOD *aod);


void runBMDSContAnalysis(struct continuous_analysis *anal, struct continuous_model_result *res, struct BMDS_results *bmdsRes, struct continuous_AOD *aod, struct continuous_GOF *gof, bool detectAdvDir);
#ifdef __cplusplus
}
#endif
