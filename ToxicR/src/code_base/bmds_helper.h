#ifdef __cplusplus
  #include <Eigen/Dense>
#endif
//#include "bmdStruct.h"
#include "dichotomous_entry_code.h"
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

const double BMDS_EPS = 1.0e-6;


//c entry
#ifdef __cplusplus
extern "C" {
#endif
//void collect_dicho_bmd_values(double *bmd_dist, struct BMD_results *BMDres);
void collect_dicho_bmd_values(struct dichotomous_analysis *anal, struct dichotomous_model_result *res, struct BMDS_results *BMDres);

void collect_cont_bmd_values(struct continuous_analysis *anal, struct continuous_model_result *res, struct BMDS_results *BMDres);


void runBMDSDichoAnalysis(struct dichotomous_analysis *anal, struct dichotomous_model_result *res, struct dichotomous_PGOF_result *gofRes, struct BMDS_results *bmdsRes);


void runBMDSContAnalysis(struct continuous_analysis *anal, struct continuous_model_result *res, struct BMDS_results *bmdsRes);
#ifdef __cplusplus
}
#endif
