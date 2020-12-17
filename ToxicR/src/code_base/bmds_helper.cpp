#include <stdio.h>
#include "bmds_helper.h"

int checkForBoundedParms(int nparms, double *parms, double *prior, struct BMDS_results *BMDSres ){
   // First find number of bounded parms
   int bounded = 0;
   for (int i=0; i<nparms; i++){
      //5*i+4 is location of min in prior array
      //5*i+5 is location of max in prior array
      if (fabs(parms[i]-prior[5*i+4]) < BMDS_EPS || fabs(parms[i]-prior[5*i+5]) < BMDS_EPS){
         bounded++;
         BMDSres->bounded[i] = false;
      }
   }
   return bounded;
}


void calcDichoAIC(struct dichotomous_analysis *anal, struct dichotomous_model_result *res, struct BMDS_results *BMDSres){

  printf("LL=%f\n",res->max);
  printf("model_df=%f\n",res->model_df);

  bool freqModel = anal->prior[0] == 0;

  // First find number of bounded parms
  int bounded = checkForBoundedParms(anal->parms, res->parms, anal->prior, BMDSres);

  double estParmCount = res->model_df - bounded;
  //if freq then model_df should be rounded to nearest whole number
  if (freqModel)
     estParmCount = round(estParmCount);
  BMDSres->AIC = 2*(res->max + estParmCount);
}

void calcContAIC(struct continuous_analysis *anal, struct continuous_model_result *res, struct BMDS_results *BMDSres){

  printf("LL=%f\n",res->max);
  printf("model_df=%f\n",res->model_df);

  bool freqModel = anal->prior[0] == 0;

  // First find number of bounded parms
   int bounded = checkForBoundedParms(anal->parms, res->parms, anal->prior, BMDSres);
  
  
  double estParmCount = res->model_df - bounded;
  //if freq then model_df should be rounded to nearest whole number
  if (freqModel)
    estParmCount = round(estParmCount);
  
    BMDSres->AIC = 2*(res->max + estParmCount);
  
  //               //    BMDSres->AIC = -9998.0;
  
}

double findQuantileVals(double *quant, double *val, int arrSize, double target){

   double retVal = -9999.0;

   for (int i=0; i < arrSize; i++){
      if (fabs(quant[i] - target) < BMDS_EPS && std::isfinite(val[i]) ){
         retVal = val[i];
      }
   }
   return retVal;
}

void collect_dicho_bmd_values(struct dichotomous_analysis *anal, struct dichotomous_model_result *res, struct BMDS_results *BMDSres){

  int distSize = res->dist_numE*2;

  double dist[distSize/2][2];
  double quant[distSize/2];
  double val[distSize/2];

  for (int i = 0; i < distSize/2; i++){
    dist[i][1] = res->bmd_dist[i];
    val[i] = res->bmd_dist[i];
  }
  for (int i = distSize/2; i < distSize; i++){
    dist[i-distSize/2][0] = res->bmd_dist[i];
    quant[i-distSize/2] = res->bmd_dist[i];
  }

  calcDichoAIC(anal, res, BMDSres);
  BMDSres->BMD = findQuantileVals(quant, val, distSize/2, 0.50);
  BMDSres->BMDL = findQuantileVals(quant, val, distSize/2, 0.05);
  BMDSres->BMDU = findQuantileVals(quant, val, distSize/2, 0.95);

}

void collect_cont_bmd_values(struct continuous_analysis *anal, struct continuous_model_result *res, struct BMDS_results *BMDSres){

  int distSize = res->dist_numE*2;
  printf("distSize: %d\n",distSize);

  double dist[distSize/2][2];
  double quant[distSize/2];
  double val[distSize/2];

  for (int i = 0; i < distSize/2; i++){
    dist[i][1] = res->bmd_dist[i];
    val[i] = res->bmd_dist[i];
  }
  for (int i = distSize/2; i < distSize; i++){
    dist[i-distSize/2][0] = res->bmd_dist[i];
    quant[i-distSize/2] = res->bmd_dist[i];
  }

  printf("dist:\n");
  for (int i = 0; i < distSize/2; i++)
  {
     printf("quant:%f, val:%f\n",dist[i][0],dist[i][1]);
     printf("quant:%f, val:%f\n",quant[i],val[i]);
  }

  calcContAIC(anal, res, BMDSres);
  BMDSres->BMD = findQuantileVals(quant, val, distSize/2, 0.50);
  BMDSres->BMDL = findQuantileVals(quant, val, distSize/2, 0.05);
  BMDSres->BMDU = findQuantileVals(quant, val, distSize/2, 0.95);

  
}
