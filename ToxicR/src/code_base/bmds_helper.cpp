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


void runBMDSDichoAnalysis(struct dichotomous_analysis *anal, struct dichotomous_model_result *res, struct dichotomous_PGOF_result *gofRes, struct BMDS_results *bmdsRes, struct dicho_AOD *bmdsAOD){

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

  //calculate model chi^2 value
  bmdsRes->chisq = 0.0;
  for (int i=0; i<gofRes->n; i++){
    bmdsRes->chisq += gofRes->residual[i]*gofRes->residual[i];
  }

  //calculate dichtomous analysis of deviance
  calc_dichoAOD(anal, res, bmdsRes, bmdsAOD);



}


void runBMDSContAnalysis(struct continuous_analysis *anal, struct continuous_model_result *res, struct BMDS_results *bmdsRes, struct continuous_AOD *aod, struct continuous_GOF *gof, bool detectAdvDir){

  if (detectAdvDir){
    determineAdvDir(anal);
  }

  estimate_sm_laplace_cont(anal, res);

  collect_cont_bmd_values(anal, res, bmdsRes);
  
  //if not suff_stat, then convert
  struct continuous_analysis GOFanal;
  //arrays are needed for conversion to suff_stat
  double doses[anal->n];
  double means[anal->n];
  double n_group[anal->n];
  double sd[anal->n];
  bool isIncreasing;
  double BMR;
  double tail_prob;
  int disttype;
  double alpha;
  int samples;
  int degree;
  int burnin;
  int parms;
  int prior_cols;

  if (anal->suff_stat){
    GOFanal = *anal;
  } else {
    //copy analysis and convert to suff_stat
    GOFanal.doses = doses;
    GOFanal.Y = means;
    GOFanal.n_group = n_group;
    GOFanal.sd = sd;
    GOFanal.isIncreasing = isIncreasing;
    GOFanal.BMR = BMR;
    GOFanal.tail_prob = tail_prob;
    GOFanal.disttype = disttype;
    GOFanal.alpha = alpha;
    GOFanal.samples = samples;
    GOFanal.degree = degree;
    GOFanal.burnin = burnin;
    GOFanal.parms = parms;
    GOFanal.prior_cols = prior_cols;
    bmdsConvertSStat(anal, &GOFanal);    
  }
 

  continuous_expected_result GOFres;
  GOFres.n = GOFanal.n;
  GOFres.expected = new double[GOFanal.n];
  GOFres.sd = new double[GOFanal.n];
  continuous_expectation(&GOFanal, res, &GOFres);


  gof->dose = GOFanal.doses;
  gof->size = GOFanal.n_group;
  gof->estMean = GOFres.expected;
  gof->obsMean = GOFanal.Y;
  gof->estSD = GOFres.sd;
  gof->obsSD = GOFanal.sd;
  for (int i=0; i<GOFanal.n; i++){
    gof->res[i] = sqrt(gof->size[i])*(gof->obsMean[i] - gof->estMean[i]) / gof->estSD[i];
  }
  if (anal->disttype == distribution::log_normal){
    for (int i=0; i<GOFanal.n; i++){
      gof->calcMean[i] = exp(log(GOFanal.Y[i]) - log(1 + pow(GOFanal.sd[i] / GOFanal.Y[i], 2.0)) / 2);
      gof->calcSD[i] = exp(sqrt(log(1.0 + pow(GOFanal.sd[i]/GOFanal.Y[i], 2.0))));
    }
  } else {
    gof->calcMean = GOFanal.Y;
    gof->calcSD = GOFanal.sd;
  }

  calc_contAOD(anal, res, bmdsRes, aod);
  
  

}


void calc_dichoAOD(struct dichotomous_analysis *DA, struct dichotomous_model_result *res, struct BMDS_results *bmdsRes, struct dicho_AOD *bmdsAOD){
  
  struct dichotomous_aod aod;
  double A1;
  int N1;
  double A2;
  int N2;
  aod.A1 = A1;
  aod.N1 = N1;
  aod.A2 = A2;
  aod.N2 = N2;

  deviance_dichotomous(DA, &aod);

  bmdsAOD->fullLL = -1*aod.A1;
  bmdsAOD->nFull = aod.N1;
  bmdsAOD->redLL = -1*aod.A2;
  bmdsAOD->nRed = aod.N2;
  bmdsAOD->fittedLL = -1*res->max;
  int bounded = 0;
  for(int i=0; i<DA->parms;i++){
    if(bmdsRes->bounded[i]){
      bounded++;
    }
  }
  bmdsAOD->nFit = DA->parms - bounded;

  double dev;
  double df;

  bmdsAOD->devFit = dev = 2*(bmdsAOD->fullLL - bmdsAOD->fittedLL);
  bmdsAOD->dfFit = df = DA->n - bmdsAOD->nFit;
  bmdsAOD->pvFit = 1.0 -gsl_cdf_chisq_P(dev, df);


  bmdsAOD->devRed = dev = 2* (bmdsAOD->fittedLL - bmdsAOD->redLL);
  bmdsAOD->dfRed = df = DA->n-1;
  bmdsAOD->pvRed = 1.0 - gsl_cdf_chisq_P(dev, df);
  
  
 
}

void calc_contAOD(struct continuous_analysis *CA, struct continuous_model_result *res, struct BMDS_results *bmdsRes, struct continuous_AOD *aod){

  continuous_deviance CD;  

  if (CA->disttype == distribution::log_normal){
    estimate_log_normal_aod(CA, &CD);
  } else {
    estimate_normal_aod(CA, &CD);
  }


  //fill aod with results and calculate tests of interest
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

  int n_rows;

  struct continuous_analysis CAnew;
  double doses[CA->n];  //allocate memory for individual size.  Will not use entire array for summary data.
  double means[CA->n];
  double n_group[CA->n];
  double sd[CA->n];
  
  CAnew.doses = doses;
  CAnew.Y = means;
  CAnew.n_group = n_group;
  CAnew.sd = sd;

  if(!CA->suff_stat){
    bmdsConvertSStat(CA, &CAnew);
    n_rows = CAnew.n;
  } else {
    //already suff stat
    n_rows = CA->n;
  }

  Eigen::MatrixXd X(n_rows,2);
  Eigen::MatrixXd W(n_rows,n_rows);
  Eigen::VectorXd Y(n_rows);
  Eigen::VectorXd beta(2);


  if(!CA->suff_stat){
    for (int i=0; i< n_rows; i++){
      X(i,0) = 1;
      X(i,1) = CAnew.doses[i];

      W(i,i) =  CAnew.n_group[i];
      Y(i) = CAnew.Y[i];
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



void bmdsConvertSStat(struct continuous_analysis *CA, struct continuous_analysis *CAss){

  //standardize the data
  int n_rows = CA->n; 
  int n_cols = CA->suff_stat?3:1;


  if(!CA->suff_stat){
    Eigen::MatrixXd Yin(n_rows,n_cols);
    Eigen::MatrixXd Xin(n_rows,1);
    Eigen::MatrixXd SSTAT, SSTAT_LN, X, Y;
    // copy the original data
    for (int i = 0; i < n_rows; i++){
      Yin(i,0) = CA->Y[i];
      Xin(i,0) = CA->doses[i];
    }
    bool canConvert = convertSStat(Yin, Xin, &SSTAT, &SSTAT_LN, &X);
    if (canConvert){
      bool isLognormal = CA->disttype == distribution::log_normal;
      if (isLognormal) {
        //Y = cleanSuffStat(SSTAT,X,false);
        Y = cleanSuffStat(SSTAT_LN,X,true,false);
      } else {
        //Y = cleanSuffStat(SSTAT_LN,X,true);
        Y = cleanSuffStat(SSTAT,X,false,false);
      }
      n_rows = X.rows();

      for(int i=0; i<n_rows; i++){
        CAss->doses[i] = X(i);
        CAss->Y[i] = Y.col(0)[i];
        CAss->n_group[i] = Y.col(1)[i];
        CAss->sd[i] = Y.col(2)[i];
      }
      CAss->n = n_rows;
       
    } else {
      std::cout<<"Error in bmdsConvertSStat"<<std::endl;
      return;
    }
  } else {
  //already suff_stat, so just copy data to arrays
    for (int i=0; i<n_rows; i++){
      CAss->doses[i] = CA->doses[i];
      CAss->Y[i] = CA->Y[i];
      CAss->n_group[i] = CA->n_group[i];
      CAss->sd[i] = CA->sd[i];
      CAss->n = CA->n;
    }

  }
  CAss->model = CA->model;
  CAss->isIncreasing = CA->isIncreasing;
  CAss->BMR = CA->BMR;
  CAss->tail_prob = CA->tail_prob;
  CAss->disttype = CA->disttype;
  CAss->alpha = CA->alpha;
  CAss->samples = CA->samples;
  CAss->degree = CA->degree;
  CAss->burnin = CA->burnin;
  CAss->parms = CA->parms;
  CAss->prior_cols = CA->prior_cols;



}
