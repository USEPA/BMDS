#include <stdio.h>
#include "bmds_helper.h"
#include "analysis_of_deviance.h"

int checkForBoundedParms(int nparms, double *parms, double *prior, struct BMDS_results *BMDSres ){
   // First find number of bounded parms
   int bounded = 0;
   for (int i=0; i<nparms; i++){
      //5*i+4 is location of min in prior array
      //5*i+5 is location of max in prior array
      if (fabs(parms[i]-prior[3*nparms+i]) < BMDS_EPS || fabs(parms[i]-prior[4*nparms+i]) < BMDS_EPS){
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

//  printf("dist:\n");
//  for (int i = 0; i < distSize/2; i++)
//  {
//     printf("quant:%f, val:%f\n",dist[i][0],dist[i][1]);
//     printf("quant:%f, val:%f\n",quant[i],val[i]);
//  }

  calcContAIC(anal, res, BMDSres);
  BMDSres->BMD = findQuantileVals(quant, val, distSize/2, 0.50);
  BMDSres->BMDL = findQuantileVals(quant, val, distSize/2, 0.05);
  BMDSres->BMDU = findQuantileVals(quant, val, distSize/2, 0.95);

  
}


void runBMDSDichoAnalysis(struct dichotomous_analysis *anal, struct dichotomous_model_result *res, struct dichotomous_PGOF_result *gofRes, struct BMDS_results *bmdsRes){

  estimate_sm_laplace_dicho(anal, res, true);

  struct dichotomous_PGOF_data gofData;
  gofData.n = anal->n;
  gofData.Y = anal->Y;
  gofData.model = anal->model;
  gofData.model_df = res->model_df;
  gofData.est_parms = res->parms;
  gofData.doses = anal->doses;
  gofData.n_group = anal->n_group;
  gofData.parms = anal->parms; 
  

  compute_dichotomous_pearson_GOF(&gofData, gofRes);

  collect_dicho_bmd_values(anal, res, bmdsRes);

}


void runBMDSContAnalysis(struct continuous_analysis *anal, struct continuous_model_result *res, struct BMDS_results *bmdsRes, struct continuous_AOD *aod, bool detectAdvDir){

  if (detectAdvDir){
    determineAdvDir(anal);
  }

  estimate_sm_laplace_cont(anal, res);

  collect_cont_bmd_values(anal, res, bmdsRes);



  //std::cout<<"calling calc_AOD"<<std::endl;
  calc_AOD(anal, res, bmdsRes, aod);
  

 // std::cout<<"AOD vals:"<<std::endl;
 // for (int i=0; i<5; i++){
 //   std::cout<<"i:"<<i<<", LL:"<<aod->LL[i]<<", nParms:"<<aod->nParms[i]<<", AIC:"<<aod->AIC[i]<<std::endl;
 // }
 // std::cout<<"additive constant:"<<aod->addConst<<std::endl;
 // std::cout<<"TOI vals:"<<std::endl;
 // for (int i=0; i<4; i++){
 //   std::cout<<"i:"<<i<<", llRatio:"<<aod->TOI->llRatio[i]<<", DF:"<<aod->TOI->DF[i]<<", pVal:"<<aod->TOI->pVal[i]<<std::endl;
 // }
  

}

void calc_AOD(struct continuous_analysis *CA, struct continuous_model_result *res, struct BMDS_results *bmdsRes, struct continuous_AOD *aod){

  continuous_deviance CD;  

  if (CA->disttype == distribution::log_normal){
    cout << "calling estimate_log_normal_aod" << std::endl;
    estimate_log_normal_aod(CA, &CD);
  } else {
    cout << "calling estimate_normal_aod" << std::endl;
    estimate_normal_aod(CA, &CD);
  }

  std::cout << "finished fit" << std::endl;

  //fill aod with results and calculate tests of interest
  std::cout << "A1:" << CD.A1 << std::endl; 
  std::cout << "A2:" << CD.A2 << std::endl; 
  std::cout << "A3:" << CD.A3 << std::endl; 
 
  aod->LL[0] = -1*CD.A1;
  aod->nParms[0] = CD.N1;
  aod->LL[1] = -1*CD.A2;
  aod->nParms[1] = CD.N2;
  aod->LL[2] = -1*CD.A3;
  aod->nParms[2] = CD.N3;
  aod->LL[3] = -1*res->max;
  int bounded = 0;
  for(int i=0; i<CA->parms;i++){
    if(bmdsRes->bounded[i]){
      bounded++;
    }
  }
  std::cout<<"CA->parms="<<CA->parms<<std::endl;
  aod->nParms[3] = CA->parms - bounded;
  //TEMP DEBUG values
  aod->LL[4] = -88.71521848;  //CD.R
  aod->nParms[4] = 2;  //CD.NR


  //add tests of interest
  //TODO:  need to check for bounded parms in A1,A2,A3,R
  for (int i=0; i<5; i++){
    aod->AIC[i] = 2*(aod->LL[i] + aod->nParms[i]);
  }  

  if (CA->suff_stat) {
    double sumN = 0;
    for(int i=0; i<CA->n; i++){
      sumN+=CA->n_group[i];
    }
    aod->addConst = -1*sumN*log(2*M_PI)/2;
  } else {
    aod->addConst = -1*CA->n*log(2*M_PI)/2;
  }

  
  
  calcTestsOfInterest(aod);
}



void calcTestsOfInterest(struct continuous_AOD *aod){

  struct testsOfInterest *TOI = aod->TOI;
  double dev;
  int df;

  //Test #1 - A2 vs Reduced - does mean and/or variance differ across dose groups
  TOI->llRatio[0] = dev = 2 * (aod->LL[1] - aod->LL[4]);
  TOI->DF[0] = df = aod->nParms[1] - aod->nParms[4];
  TOI->pVal[0] = (dev < 0.0) ? -1 : 1.0 - gsl_cdf_chisq_P(dev, df);

  //Test #2 - A1 vs A2 - homogeneity of variance across dose groups
  TOI->llRatio[1] = dev = 2 * (aod->LL[1] - aod->LL[0]);
  TOI->DF[1] = df = aod->nParms[1] - aod->nParms[0];
  TOI->pVal[1] = (dev < 0.0) ? -1 : 1.0 - gsl_cdf_chisq_P(dev, df);

  //Test #3 - A2 vs A3 - Does the model describe variances adequately
  TOI->llRatio[2] = dev = 2 * (aod->LL[1] - aod->LL[2]);
  TOI->DF[2] = df = aod->nParms[1] - aod->nParms[2];
  TOI->pVal[2] = (dev < 0.0) ? -1 : 1.0 - gsl_cdf_chisq_P(dev, df);

  //Test #4 - A3 vs Fitted - does the fitted model describe the obs data adequately 
  TOI->llRatio[3] = dev = 2 * (aod->LL[2] - aod->LL[3]);
  TOI->DF[3] = df = aod->nParms[2] - aod->nParms[3];
  TOI->pVal[3] = (dev < 0.0) ? -1 : 1.0 - gsl_cdf_chisq_P(dev, df);
 
}

void determineAdvDir(struct continuous_analysis *CA){

  //standardize the data
  int n_rows = CA->n; int n_cols = CA->suff_stat?3:1; 
  
  Eigen::MatrixXd Yin(n_rows,n_cols); 
  Eigen::MatrixXd Xin(n_rows,1); 
  Eigen::MatrixXd SSTAT, SSTAT_LN, UX;
  Eigen::MatrixXd Y_LN, Y_N;

  if(!CA->suff_stat){
  // copy the original data
    for (int i = 0; i < n_rows; i++){
      Yin(i,0) = CA->Y[i]; 
      Xin(i,0) = CA->doses[i]; 
      if(CA->suff_stat){
        
        Yin(i,2) = CA->sd[i]; 
        Yin(i,1) = CA->n_group[i]; 
      }
    }

    bool canConvert = convertSStat(Yin, Xin, &SSTAT, &SSTAT_LN, &UX);
    if (canConvert){
      Y_N = cleanSuffStat(SSTAT,UX,false);  
      Y_LN = cleanSuffStat(SSTAT_LN,UX,true); 
      n_rows = UX.rows();

    } else {
      std::cout<<"Error in determineAdvDir"<<std::endl;
      return;
    }
  }



  Eigen::MatrixXd X(n_rows,2);
  Eigen::MatrixXd W(n_rows,n_rows);
  Eigen::VectorXd Y(n_rows);
  Eigen::VectorXd beta(2);


  if(!CA->suff_stat){
   
      for (int i=0; i< n_rows; i++){
        X(i,0) = 1;
        X(i,1) = UX(i);
  
        W(i,i) =  Y_N(i,1);
        Y(i) = Y_N(i,0);
    }  


  } else {
    for (int i=0; i< n_rows; i++){
      X(i,0) = 1;
      X(i,1) = CA->doses[i];

      W(i,i) =  CA->n_group[i];
      Y(i) = CA->Y[i];
    }
 
  }

  beta = (X.transpose() * W * X).colPivHouseholderQr().solve(X.transpose()*W*Y);

  if (beta(1) > 0 ) {
    CA->isIncreasing = true;
  } else {
    CA->isIncreasing = false;
  }

}

