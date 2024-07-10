#ifdef WIN32
	#include "pch.h"
#else
	#include "stdafx.h"
#endif
#include <stdio.h>
#include <math.h>
#include <ctime>
//#include <cmath>
#include "bmds_helper.h"
#include "analysis_of_deviance.h"

// calendar versioning; see https://peps.python.org/pep-0440/#pre-releases
std::string BMDS_VERSION = "2023.10a1";

int checkForBoundedParms(int nparms, double *parms, double *lowerBound, double *upperBound, struct BMDS_results *BMDSres ){
   // First find number of bounded parms
   int bounded = 0;
   BMDSres->bounded.clear();
   for (int i=0; i<nparms; i++){
      BMDSres->bounded.push_back(false);
      //5*i+4 is location of min in prior array
      //5*i+5 is location of max in prior array 
      if (fabs(parms[i]-lowerBound[i]) < BMDS_EPS || fabs(parms[i]-upperBound[i]) < BMDS_EPS){
         bounded++;
         BMDSres->bounded[i] = true;
      }
   }
   return bounded;
}


void calcDichoAIC(struct dichotomous_analysis *anal, struct dichotomous_model_result *res, struct BMDS_results *BMDSres, double estParmCount){

  bool freqModel = anal->prior[0] == 0;

//  // First find number of bounded parms
  //int bounded = checkForBoundedParms(anal->parms, res->parms, anal->prior, BMDSres);

  double* lowerBound = (double*)malloc(anal->parms * sizeof(double));
  double* upperBound = (double*)malloc(anal->parms * sizeof(double));

  //copy prior values
  for (int i=0; i<anal->parms; i++){
    lowerBound[i] = anal->prior[3*anal->parms+i];
    upperBound[i] = anal->prior[4*anal->parms+i];
  }

  //scale priors as needed
  rescale_dichoParms(anal->model, lowerBound);
  rescale_dichoParms(anal->model, upperBound);

  int bounded = checkForBoundedParms(anal->parms, res->parms, lowerBound, upperBound, BMDSres);

  estParmCount = anal->parms - bounded;

  //if freq then model_df should be rounded to nearest whole number
  if (freqModel)
     estParmCount = round(estParmCount);
  BMDSres->AIC = 2*(res->max + estParmCount);

  free(lowerBound);
  free(upperBound);
}

void calcContAIC(struct continuous_analysis *anal, struct continuous_model_result *res, struct BMDS_results *BMDSres){

  bool freqModel = anal->prior[0] == 0;

  // First find number of bounded parms
  double* lowerBound = (double*)malloc(anal->parms * sizeof(double));
  double* upperBound = (double*)malloc(anal->parms * sizeof(double));

  //copy prior values
  for (int i=0; i<anal->parms; i++){
     lowerBound[i] = anal->prior[3*anal->parms+i];
     upperBound[i] = anal->prior[4*anal->parms+i];
  }
  //scale priors as needed
  rescale_contParms(anal, lowerBound);
  rescale_contParms(anal, upperBound); 
  
  if (anal->model == cont_model::exp_3){
    for (int i=2; i<anal->parms-1; i++){
      lowerBound[i] = lowerBound[i+1];
      upperBound[i] = upperBound[i+1];
    } 
  }
  int bounded = checkForBoundedParms(res->nparms, res->parms, lowerBound, upperBound, BMDSres);
  //int bounded = checkForBoundedParms(anal->parms, res->parms, lowerBound, upperBound, BMDSres);
  
  double estParmCount = res->model_df - bounded;
  //if freq then model_df should be rounded to nearest whole number
  if (freqModel)
    estParmCount = round(estParmCount);
  
    BMDSres->AIC = 2*(res->max + estParmCount);
  
  //               //    BMDSres->AIC = -9998.0;
  free(lowerBound);
  free(upperBound);
  
}

double findQuantileVals(double *quant, double *val, int arrSize, double target){

   double retVal = BMDS_MISSING;

   for (int i=0; i < arrSize; i++){
      if (fabs(quant[i] - target) < BMDS_EPS && std::isfinite(val[i]) ){
         //exact match
         retVal = val[i];
         break;
      } else if (quant[i] > target && i>0){ 
        // linear interpolation
        retVal = val[i-1] + ((val[i] - val[i-1])/(quant[i] - quant[i-1])) * (target - quant[i-1]);
        break;
      }
   }
   return retVal;
}

void collect_dicho_bmd_values(struct dichotomous_analysis *anal, struct dichotomous_model_result *res, struct BMDS_results *BMDSres, double estParmCount){

  int distSize = res->dist_numE*2;

  double* quant = (double*)malloc(distSize / 2 * sizeof(double));
  double* val = (double*)malloc(distSize / 2 * sizeof(double));

  for (int i = 0; i < distSize/2; i++){
    val[i] = res->bmd_dist[i];
  }
  for (int i = distSize/2; i < distSize; i++){
    quant[i-distSize/2] = res->bmd_dist[i];
  }

  calcDichoAIC(anal, res, BMDSres, estParmCount);
  BMDSres->BMD = findQuantileVals(quant, val, distSize/2, 0.50);
  BMDSres->BMDL = findQuantileVals(quant, val, distSize/2, 0.05);
  BMDSres->BMDU = findQuantileVals(quant, val, distSize/2, 0.95);

  free(quant);
  free(val);

}


void collect_dichoMA_bmd_values(struct dichotomousMA_analysis *anal, struct dichotomousMA_result *res, struct BMDSMA_results *BMDSres){

  int distSize = res->dist_numE*2;

  double* quantMA = (double*)malloc(distSize / 2 * sizeof(double));
  double* valMA = (double*)malloc(distSize / 2 * sizeof(double));

  for (int i = 0; i < distSize/2; i++){
    valMA[i] = res->bmd_dist[i];
  }
  for (int i = distSize/2; i < distSize; i++){
    quantMA[i-distSize/2] = res->bmd_dist[i];
  }

//  calculate MA quantiles
  BMDSres->BMD_MA = findQuantileVals(quantMA, valMA, distSize/2, 0.50);
  BMDSres->BMDL_MA = findQuantileVals(quantMA, valMA, distSize/2, 0.05);
  BMDSres->BMDU_MA = findQuantileVals(quantMA, valMA, distSize/2, 0.95);

// calculate individual model quantiles
  for (int j=0; j<anal->nmodels; j++){
      for (int i = 0; i < distSize/2; i++){
        valMA[i] = res->models[j]->bmd_dist[i];
      }
      for (int i = distSize/2; i < distSize; i++){
        quantMA[i-distSize/2] = res->models[j]->bmd_dist[i];
      }
      BMDSres->BMD[j] = findQuantileVals(quantMA, valMA, distSize/2, 0.50);
      BMDSres->BMDL[j] = findQuantileVals(quantMA, valMA, distSize/2, 0.05);
      BMDSres->BMDU[j] = findQuantileVals(quantMA, valMA, distSize/2, 0.95);
  }
  free(quantMA);
  free(valMA);
}

void collect_cont_bmd_values(struct continuous_analysis *anal, struct continuous_model_result *res, struct BMDS_results *BMDSres){

  int distSize = res->dist_numE*2;

  double* contVal = (double*)malloc(distSize / 2 * sizeof(double));
  double* contQuant = (double*)malloc(distSize / 2 * sizeof(double));

  for (int i = 0; i < distSize/2; i++){
    contVal[i] = res->bmd_dist[i];
  }

  for (int i = distSize/2; i < distSize; i++){
    contQuant[i-distSize/2] = res->bmd_dist[i];
  }

  calcContAIC(anal, res, BMDSres);
  BMDSres->BMD = findQuantileVals(contQuant, contVal, distSize/2, 0.50);
  BMDSres->BMD = findQuantileVals(contQuant, contVal, distSize/2, 0.50);
  BMDSres->BMDL = findQuantileVals(contQuant, contVal, distSize/2, 0.05);
  BMDSres->BMDU = findQuantileVals(contQuant, contVal, distSize/2, 0.95);

  free(contVal);
  free(contQuant);
}


void BMDS_ENTRY_API __stdcall runBMDSDichoAnalysis(struct dichotomous_analysis *anal, struct dichotomous_model_result *res, struct dichotomous_GOF *gof, struct BMDS_results *bmdsRes, struct dicho_AOD *bmdsAOD){
 
//  std::cout<<"degree = "<<anal->degree;
//  std::cout<<"degree:"<<anal->degree<<std::endl;
//        for (int k=0; k<anal->prior_cols*anal->parms; k++){
//          std::cout<<anal->prior[k]<<", ";
//        }
//        std::cout<<std::endl;
 
  bmdsRes->validResult = false;
  bmdsRes->bounded.resize(anal->parms);
  fill(bmdsRes->bounded.begin(), bmdsRes->bounded.end(), false);
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

  struct dichotomous_PGOF_result gofRes;
  double* gofExpected = (double*)malloc(anal->n * sizeof(double));
  double* gofResidual = (double*)malloc(anal->n * sizeof(double));
  double gofTestStat = BMDS_MISSING;
  double gofPVal = BMDS_MISSING;
  double gofDF = BMDS_MISSING;
  double* ebUpper = (double*)malloc(anal->n * sizeof(double));
  double* ebLower = (double*)malloc(anal->n * sizeof(double));
  gofRes.n = anal->n;
  gofRes.expected = gofExpected;
  gofRes.residual = gofResidual;
  gofRes.test_statistic = gofTestStat;

  gofRes.p_value = gofPVal;
  gofRes.df = gofDF;
  
  compute_dichotomous_pearson_GOF(&gofData, &gofRes);
  
  gof->test_statistic = gofRes.test_statistic;

  //these will be updated later with bounded parm info
  gof->p_value = gofRes.p_value;
  gof->df = gofRes.df;


  gof->n = gofRes.n;
  for (int i=0; i<gofRes.n; i++){
    gof->expected.push_back(gofRes.expected[i]);
    gof->residual.push_back(gofRes.residual[i]);
  }
  
  
  //do error bar calcs
  double pHat;
  double z;
  double eb1;
  double eb2;
  double ebDenom;
  double gofAlpha = 0.05;  //Alpha value for 95% confidence limit
  z = gsl_cdf_ugaussian_Pinv(1.0 - gofAlpha / 2); //Z score
  for (int i=0; i<gof->n; i++){
    pHat = anal->Y[i]/anal->n_group[i];  //observed probability
    eb1 = (2 * anal->Y[i] + z*z);
    eb2 = z*sqrt(z*z - (2 + 1/anal->n_group[i]) + 4*pHat*((anal->n_group[i] - anal->Y[i]) + 1));
    ebDenom = 2.0 * (anal->n_group[i] + z * z);
    gof->ebLower.push_back((eb1 - 1 - eb2)/ ebDenom); 
    gof->ebUpper.push_back((eb1 + 1 + eb2)/ ebDenom);
  }

  //calculate model chi^2 value
  bmdsRes->chisq = 0.0;
  for (int i=0; i<gofRes.n; i++){
    bmdsRes->chisq += gofRes.residual[i]*gofRes.residual[i];
  }

  //calculate bayesian BIC_equiv
  Eigen::MatrixXd cov(res->nparms,res->nparms);
  int row = 0;
  int col = 0;
  for(int i=0; i<res->nparms*res->nparms; i++){
    col = i/res->nparms;
    row = i - col*res->nparms;
    cov(row,col) = res->cov[i];
  }

  bmdsRes->BIC_equiv = res->nparms / 2.0 *log(2.0*M_PI) + res->max + 0.5*log(max(0.0, cov.determinant()));
  bmdsRes->BIC_equiv = -1*bmdsRes->BIC_equiv;


  //calculate dichtomous analysis of deviance
  struct dichotomous_aod aod;
  double A1 = BMDS_MISSING;
  int N1 = BMDS_MISSING;
  double A2 = BMDS_MISSING;
  int N2 = BMDS_MISSING;
  aod.A1 = A1;
  aod.N1 = N1;
  aod.A2 = A2;
  aod.N2 = N2;
  
  calc_dichoAOD(anal, res, bmdsRes, bmdsAOD, &aod);

  rescale_dichoParms(anal->model, res->parms);

  double estParmCount = 0;
  collect_dicho_bmd_values(anal, res, bmdsRes, estParmCount);

  //incorporate affect of bounded parameters
  int bounded = 0;
  for (int i=0; i<anal->parms; i++){
    if (bmdsRes->bounded[i]){
      bounded++;
    }
  }

  bmdsAOD->nFit = anal->parms - bounded;  //number of estimated parameter
  bmdsAOD->dfFit = anal->n - bmdsAOD->nFit;  //nObs - nEstParms 

  if (bmdsAOD->devFit < 0 || bmdsAOD->dfFit < 0) {
    bmdsAOD->pvFit = BMDS_MISSING;
  } else {
    bmdsAOD->pvFit = 1.0 -gsl_cdf_chisq_P(bmdsAOD->devFit, bmdsAOD->dfFit);
  }


  //update df for frequentist models only
  if(anal->prior[1] == 0){
    //frequentist
    gofRes.df = bmdsAOD->dfFit;
  } 
  gof->df = gofRes.df;

  if ( gof->df > 0.0){
    gofRes.p_value        =   1.0 - gsl_cdf_chisq_P(bmdsRes->chisq, gof->df);
  }else{
    gofRes.p_value       = 1.0;
  }

  //gofRes.p_value = 
  if (gof->df <= 0.0){
    gof->p_value = BMDS_MISSING;
  } else {
    gof->p_value = gofRes.p_value;
  }


  for (int i=0; i< anal->parms; i++){
    //std err is sqrt of covariance diagonals unless parameter hit a bound, then report NA
    bmdsRes->stdErr.push_back(BMDS_MISSING);
    bmdsRes->lowerConf.push_back(BMDS_MISSING);
    bmdsRes->upperConf.push_back(BMDS_MISSING);
  } 

  calcParmCIs_dicho (res, bmdsRes);


  //compare Matt's BMD to BMD from CDF
  if (abs(bmdsRes->BMD - res->bmd) > BMDS_EPS){
    std::cout<<"Warning: CDF BMD differs from model BMD" << std::endl;
  }
  //Use Matt's BMD by default
  bmdsRes->BMD = res->bmd;

  clean_dicho_results(res, gof, bmdsRes, bmdsAOD);

  bool goodCDF = false;
  for (int i=0; i<res->dist_numE;i++){
    if (res->bmd_dist[i] != BMDS_MISSING){
      goodCDF = true;
    }
  }

  if (bmdsRes->BMD != BMDS_MISSING && !std::isinf(bmdsRes->BMD) && std::isfinite(bmdsRes->BMD) && goodCDF) {
    bmdsRes->validResult = true;
  }

}

void BMDS_ENTRY_API __stdcall runBMDSDichoMA(struct dichotomousMA_analysis *MA, struct dichotomous_analysis *DA,  struct dichotomousMA_result *res, struct BMDSMA_results *bmdsRes){

  estimate_ma_laplace_dicho(MA, DA, res);


  collect_dichoMA_bmd_values(MA, res, bmdsRes);
  for(int i=0; i<MA->nmodels; i++){
    rescale_dichoParms(MA->models[i], res->models[i]->parms);
  }

  //calc error bars for graphs
  double pHat;
  double z;
  double eb1;
  double eb2;
  double ebDenom;
  double gofAlpha = 0.05;  //Alpha value for 95% confidence limit
  z = gsl_cdf_ugaussian_Pinv(1.0 - gofAlpha / 2); //Z score
    for (int i=0; i<DA->n; i++){
    pHat = DA->Y[i]/DA->n_group[i];  //observed probability
    eb1 = (2 * DA->Y[i] + z*z);
    eb2 = z*sqrt(z*z - (2 + 1/DA->n_group[i]) + 4*pHat*((DA->n_group[i] - DA->Y[i]) + 1));
    ebDenom = 2.0 * (DA->n_group[i] + z * z);
    bmdsRes->ebLower[i] = (eb1 - 1 - eb2)/ ebDenom;
    bmdsRes->ebUpper[i] = (eb1 + 1 + eb2)/ ebDenom;
  }

  clean_dicho_MA_results(res, bmdsRes);

}

//transform parameters or priors which were computed using a logistic dist.
void rescale_dichoParms(int model, double *parms){
  //rescale background parameters for all models and v parameter for DHill model
  switch(model){
    case dich_model::d_multistage:
    case dich_model::d_weibull:
    case dich_model::d_gamma:
    case dich_model::d_loglogistic:
    case dich_model::d_qlinear:
    case dich_model::d_logprobit:
      //rescale background parameter
      parms[0] = 1.0/(1.0+exp(-1.0*parms[0]));
      break;
    case dich_model::d_hill:
      //rescale background and v parameter
      parms[0] = 1.0/(1.0+exp(-1.0*parms[0]));
      parms[1] = 1.0/(1.0+exp(-1.0*parms[1]));
      break;
    default:
      break;
  }

}

void calcParmCIs_dicho (struct dichotomous_model_result *res, struct BMDS_results *bmdsRes) {

  double *adj = new double[res->nparms];  //adjustments for transformed parameters
  int diagInd = 0;  //1d index of 2d diagonal location
  double tmp = 0;  //used for tmp calculations

  

  for (int i=0; i<res->nparms; i++){
    diagInd = i*(1+res->nparms); //gives index of diagonal for 1d flattened array
    adj[i] = 1.0;
    if (!bmdsRes->bounded[i]){
      bmdsRes->stdErr[i] = sqrt(res->cov[diagInd]); 
    }
  }      

  //now check if parameters were transformed and modify adj values accordingly
  //if transformed, then multiply adj * derivative of transform (at parm value)
  switch(res->model){
    //transform parameters which were computed using a logistic dist.
    //void rescale_dichoParms(int model, struct dichotomous_model_result *res){
    case dich_model::d_multistage:
    case dich_model::d_weibull:
    case dich_model::d_gamma:
    case dich_model::d_loglogistic:
    case dich_model::d_qlinear:
    case dich_model::d_logprobit:
      //background parameter
      tmp = exp(-1*res->parms[0]);
      tmp = tmp/((1.0+tmp)*(1.0+tmp));
      adj[0] = tmp*tmp;
      break;
    case dich_model::d_hill:
      //background and v parameter
      tmp = exp(-1*res->parms[0]);     
      tmp = tmp/((1.0+tmp)*(1.0+tmp));
      adj[0] = tmp*tmp;
      tmp = exp(-1*res->parms[1]);
      tmp = tmp/((1.0+tmp)*(1.0+tmp));
      adj[1] = tmp*tmp;
      break;
    default:
      break;
  }

  //now calculate upper and lower confidence intervals
   
  for (int i=0; i<res->nparms; i++){
    if (bmdsRes->bounded[i]){
      bmdsRes->lowerConf[i] = BMDS_MISSING;
      bmdsRes->upperConf[i] = BMDS_MISSING;
    } else {
      bmdsRes->stdErr[i] = bmdsRes->stdErr[i] * adj[i];  
      bmdsRes->lowerConf[i] = res->parms[i] - bmdsRes->stdErr[i]*BMDS_QNORM;
      bmdsRes->upperConf[i] = res->parms[i] + bmdsRes->stdErr[i]*BMDS_QNORM;
    }
  }
 
  delete [] adj; 
 
}


void BMDS_ENTRY_API __stdcall runBMDSContAnalysis(struct continuous_analysis *anal, struct continuous_model_result *res, struct BMDS_results *bmdsRes, struct continuous_AOD *aod, struct continuous_GOF *gof, bool *detectAdvDir, bool *restricted){
  bmdsRes->validResult = false;
  anal->transform_dose = false;
  //if (anal->model == cont_model::polynomial && anal->disttype == distribution::log_normal){
  if (anal->model != cont_model::exp_3 && anal->model != cont_model::exp_5){
    if(anal->disttype == distribution::log_normal){
      return; 
    }
  }
  if (*detectAdvDir){
    determineAdvDir(anal);
    int ind;
    if (*restricted) {
      switch(anal->model){
        case cont_model::exp_3:
        case cont_model::exp_5:
          if (anal->prior[0] == 0){   //checks if frequentist model
            if(anal->isIncreasing){
              if (anal->disttype == distribution::normal_ncv){
                anal->prior[20] = 0.0;  //c min
              } else {
                anal->prior[17] = 0.0;  //c min
              }
            } else {
              if (anal->disttype == distribution::normal_ncv){
                anal->prior[26] = 0.0;  //c max
              } else {
                anal->prior[22] = 0.0;  //c max
              }
            }
          }
          break;
        case cont_model::polynomial:
          if(anal->prior[0] == 0){  //checks if frequentist model
            int numRows = 2+anal->degree;
            int ind; //index in prior array
            if (anal->disttype == distribution::normal_ncv){
              numRows++;
            }
            //set ind to target min or max for 1st beta parameter then adjust as needed
            //min for increasing, max for decreasing
            if(anal->isIncreasing){
              ind = 3*numRows+1;          
            } else { 
              ind = 4*numRows+1;
            }  
            //beta terms
            for (int i=ind; i<ind+anal->degree; i++){
              anal->prior[i] = 0.0;  // beta min or max 
            }
          }
          break;
      }  //end switch
    } //end if restricted
  } //end if detectAdvDir

  //need to handle issue with decreasing response datasets for relative deviation
  //easiest workaround is to modify BMRF before sending to estimate_sm_laplace
  if (anal->BMD_type == 3) {  //rel dev
    if(!anal->isIncreasing){
      anal->BMR = 1.0 - anal->BMR;
    }
  }

  estimate_sm_laplace_cont(anal, res);

  //if not suff_stat, then convert
  struct continuous_analysis GOFanal;
  //arrays are needed for conversion to suff_stat
  double* doses = (double*)malloc(anal->n*sizeof(double));
  double* means = (double*)malloc(anal->n * sizeof(double));
  double* n_group = (double*)malloc(anal->n * sizeof(double));
  double* sd = (double*)malloc(anal->n * sizeof(double));
  for (int i=0; i<anal->n; i++){
     means[i] = anal->Y[i];
     doses[i] = anal->doses[i];
  }
  bool isIncreasing = true;
  double BMR = BMDS_MISSING;
  double tail_prob = BMDS_MISSING;
  int disttype = BMDS_MISSING;
  double alpha = BMDS_MISSING;
  int samples = BMDS_MISSING;
  int degree = BMDS_MISSING;
  int burnin = BMDS_MISSING;
  int parms = BMDS_MISSING;
  int prior_cols = BMDS_MISSING;

  if (anal->suff_stat){
    GOFanal = *anal;
  } else {
    //copy analysis and convert to suff_stat
    GOFanal.n = anal->n;
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
    GOFanal.disttype = distribution::normal;  //needed for all distrubutions to avoid error in ln conversion to suff_stats
    bmdsConvertSStat(anal, &GOFanal, true);    
  }
 

  continuous_expected_result GOFres;
  GOFres.n = GOFanal.n;
  GOFres.expected = new double[GOFanal.n];
  GOFres.sd = new double[GOFanal.n];

  continuous_expectation(&GOFanal, res, &GOFres);

  for (int i=0; i<GOFanal.n; i++){
    gof->dose.push_back(GOFanal.doses[i]);
    gof->size.push_back(GOFanal.n_group[i]);
    gof->estMean.push_back(GOFres.expected[i]);
    gof->obsMean.push_back(GOFanal.Y[i]);
    gof->estSD.push_back(GOFres.sd[i]);
    gof->obsSD.push_back(GOFanal.sd[i]);
    gof->res.push_back(sqrt(gof->size[i])*(gof->obsMean[i] - gof->estMean[i]) / gof->estSD[i]);
    gof->n = GOFanal.n;
  }
  if (anal->disttype == distribution::log_normal){
    for (int i=0; i<GOFanal.n; i++){
      gof->calcMean.push_back(exp(log(GOFanal.Y[i]) - log(1 + pow(GOFanal.sd[i] / GOFanal.Y[i], 2.0)) / 2));
      gof->calcSD.push_back(exp(sqrt(log(1.0 + pow(GOFanal.sd[i]/GOFanal.Y[i], 2.0)))));
    }
    for (int i=0; i<GOFanal.n; i++){
      gof->estMean.push_back(exp(gof->estMean[i]+pow(exp(res->parms[res->nparms-1]),2)/2));
      gof->res.push_back(sqrt(gof->size[i])*(gof->obsMean[i] - gof->estMean[i]) / gof->estSD[i]);
    }
  } else {
	for (int i=0; i<GOFanal.n; i++){
	  gof->calcMean.push_back(GOFanal.Y[i]);
          gof->calcSD.push_back(GOFanal.sd[i]);
	}
  }

  double ebUpper, ebLower;
  for (int i=0; i<GOFanal.n; i++){
    ebLower = gof->calcMean[i] + gsl_cdf_tdist_Pinv(0.025, gof->n - 1) * (gof->obsSD[i]/sqrt(gof->size[i]));
    ebUpper = gof->calcMean[i] + gsl_cdf_tdist_Pinv(0.975, gof->n - 1) * (gof->obsSD[i]/sqrt(gof->size[i]));
    gof->ebLower.push_back(ebLower);
    gof->ebUpper.push_back(ebUpper);
  }
 //calculate bayesian BIC_equiv
  Eigen::MatrixXd cov(res->nparms,res->nparms);
  int row = 0;
  int col = 0;
  for(int i=0; i<res->nparms*res->nparms; i++){
    col = i/res->nparms;
    row = i - col*res->nparms;
    cov(row,col) = res->cov[i];
  }
  
  bmdsRes->BIC_equiv = res->nparms / 2.0 *log(2.0*M_PI) + res->max + 0.5*log(max(0.0, cov.determinant()));
  bmdsRes->BIC_equiv = -1*bmdsRes->BIC_equiv;
  collect_cont_bmd_values(anal, res, bmdsRes);
  
  aod->LL.resize(5);
  aod->nParms.resize(5);
  aod->AIC.resize(5);  
  aod->TOI.llRatio.resize(4);
  aod->TOI.DF.resize(4);
  aod->TOI.pVal.resize(4);
  calc_contAOD(anal, &GOFanal, res, bmdsRes, aod);

  rescale_contParms(anal, res->parms); 

  for (int i=0; i< anal->parms; i++){
    bmdsRes->stdErr.push_back(BMDS_MISSING);
    bmdsRes->lowerConf.push_back(BMDS_MISSING);
    bmdsRes->upperConf.push_back(BMDS_MISSING);
  }

  calcParmCIs_cont(res, bmdsRes);

  //compare Matt's BMD to BMD from CDF
  if (abs(bmdsRes->BMD - res->bmd) > BMDS_EPS){
      std::cout<<"Warning: CDF BMD differs from model BMD" << std::endl;
  }
  //Use Matt's BMD by default
  bmdsRes->BMD = res->bmd;

  clean_cont_results(res, bmdsRes, aod, gof);

  bool goodCDF = false;
  for (int i=0; i<res->dist_numE;i++){
    if (res->bmd_dist[i] != BMDS_MISSING){
      goodCDF = true;
    }
  }   

  if (bmdsRes->BMD != BMDS_MISSING && !std::isinf(bmdsRes->BMD) && std::isfinite(bmdsRes->BMD) && goodCDF) { 
    bmdsRes->validResult = true;
  }
 // //compute confidence number
 // std::cout<<"cov: " << cov << std::endl;
 // Eigen::VectorXcd eivals = cov.eigenvalues();
 // std::cout << "Eigenvalues are: " << std::endl << eivals << std::endl;
 // Eigen::VectorXd eVals = eivals.real();
 // double min = *std::min_element(std::begin(eVals.array()), std::end(eVals.array()));
 // double max = *std::max_element(std::begin(eVals.array()), std::end(eVals.array()));
 // std::cout << "Min: " << min << std::endl;
 // std::cout << "Max: " << max << std::endl;
 // std::cout << "Condition number: " << max/min << std::endl;
 delete [] GOFres.expected;
 delete [] GOFres.sd;
}


void rescale_contParms(struct continuous_analysis *CA, double *parms){
  //rescale alpha parameter for all continuous models
  //assumes alpha parameter is always last
  switch (CA->model){
    case cont_model::hill:
    case cont_model::power:
    case cont_model::funl:
    case cont_model::polynomial:
      //rescale alpha
      parms[CA->parms-1] = exp(parms[CA->parms-1]); 
      break;
//    case cont_model::exp_3:
//      //rescale g for log_normal
//      if (CA->disttype == distribution::log_normal){
//        parms[0] = parms[0]*exp(pow(exp(parms[CA->parms-2]),2)/2);
//      }
//      break;
    case cont_model::exp_5:
//      //rescale g for log_normal
//      if (CA->disttype == distribution::log_normal){
//        parms[0] = parms[0]*exp(pow(exp(parms[CA->parms-2]),2)/2);
//      }
      //rescale c
      parms[2] = exp(parms[2]);
      break;
  }
}

void calcParmCIs_cont(struct continuous_model_result *res, struct BMDS_results *bmdsRes){

  double *adj = new double[res->nparms];  //adjustments for transformed parameters
  int diagInd = 0;  //1d index of 2d diagonal location
  double tmp = 0;  //used for tmp calculations

  for (int i=0; i<res->nparms; i++){
    diagInd = i*(1+res->nparms); //gives index of diagonal for 1d flattened array
    adj[i] = 1.0;
    if (!bmdsRes->bounded[i]){
      bmdsRes->stdErr[i] = sqrt(res->cov[diagInd]);
    }
  }

  //now check if parameters were transformed and modify adj values accordingly
  //if transformed, then multiply adj * derivative of transform (at parm value)
  switch(res->model){  
    case cont_model::hill:
    case cont_model::power:
    case cont_model::funl:
    case cont_model::polynomial:
      //rescale alpha  //at this point alpha has already been rescaled, so no need to use exp(parm)
      //since actual parm is ln(alpha), then gprime(theta) = exp(theta) = exp(ln(alpha)) = alpha
      //tmp = exp(res->parms[res->nparms-1]);
      tmp = res->parms[res->nparms-1];
      adj[res->nparms-1] = adj[res->nparms-1]*tmp*tmp;
      break;
    case cont_model::exp_5:
      //rescale c
      //tmp = exp(res->parms[2]); 
      //adj[2] = adj[2]*tmp*tmp;
      break;
    default:
      break;  
  }

  //now calculate upper and lower confidence intervals
  for (int i=0; i<res->nparms; i++){
    if (bmdsRes->bounded[i]){
      bmdsRes->lowerConf[i] = BMDS_MISSING;
      bmdsRes->upperConf[i] = BMDS_MISSING;
    } else {
      bmdsRes->stdErr[i] = bmdsRes->stdErr[i] * adj[i];
      bmdsRes->lowerConf[i] = res->parms[i] - bmdsRes->stdErr[i]*BMDS_QNORM;
      bmdsRes->upperConf[i] = res->parms[i] + bmdsRes->stdErr[i]*BMDS_QNORM;
    }
  }

}


void calc_dichoAOD(struct dichotomous_analysis *DA, struct dichotomous_model_result *res, struct BMDS_results *bmdsRes, struct dicho_AOD *bmdsAOD, struct dichotomous_aod *aod){
  
  deviance_dichotomous(DA, aod);

  bmdsAOD->fullLL = -1*aod->A1;
  bmdsAOD->nFull = aod->N1;
  bmdsAOD->redLL = -1*aod->A2;
  bmdsAOD->nRed = aod->N2;
  bmdsAOD->fittedLL = -1*res->max;
  int bounded = 0;
  for(int i=0; i<DA->parms;i++){
    if(bmdsRes->bounded[i]){
      bounded++;
    }
  }

  bmdsAOD->devFit = 2*(bmdsAOD->fullLL - bmdsAOD->fittedLL);

  //these fitted model values will be calculated later after bounded parms have been determined
  bmdsAOD->nFit = BMDS_MISSING;
  bmdsAOD->dfFit = BMDS_MISSING;
  bmdsAOD->pvFit = BMDS_MISSING;

  double dev;
  double df;

  bmdsAOD->devRed = dev = 2* (bmdsAOD->fullLL - bmdsAOD->redLL);
  bmdsAOD->dfRed = df = DA->n-1;
  if (dev < 0 || df < 0) {
    bmdsAOD->pvRed = BMDS_MISSING;
  } else {
    bmdsAOD->pvRed = 1.0 - gsl_cdf_chisq_P(dev, df);
  }

}

void calc_contAOD(struct continuous_analysis *CA, struct continuous_analysis *GOFanal, struct continuous_model_result *res, struct BMDS_results *bmdsRes, struct continuous_AOD *aod){
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
  //A3 is strictly NCV test.  If CV use results from A1
  if (CA->disttype == distribution::normal_ncv) {
    aod->LL[2] = -1*CD.A3;
    aod->nParms[2] = CD.N3;
  } else {
    aod->LL[2] = -1*CD.A1;
    aod->nParms[2] = CD.N1;
  }
  aod->LL[3] = -1*res->max;
  int bounded = 0;
  for(int i=0; i<CA->parms;i++){
    if(bmdsRes->bounded[i]){
      bounded++;
    }
  }
  //aod->nParms[3] = CA->parms - bounded;
  aod->nParms[3] = res->nparms - bounded;
  aod->LL[4] = -1*CD.R;
  aod->nParms[4] = CD.NR;

  //add tests of interest
  //TODO:  need to check for bounded parms in A1,A2,A3,R
  double sumN = 0;
  for (int i=0; i<5; i++){
    aod->AIC[i] = 2*(-1*aod->LL[i] + aod->nParms[i]);
  }  

  if (CA->suff_stat) {
    for(int i=0; i<CA->n; i++){
      sumN+=CA->n_group[i];
    }
  } else {
    sumN = CA->n;
  }
  aod->addConst = -1*sumN*log(2*M_PI)/2;
  if (CA->disttype == distribution::log_normal) {
    double tmp = 0;
    if (CA->suff_stat){
       for (int i=0; i<CA->n;i++){
          tmp += (log(CA->Y[i])-log(1+pow((CA->sd[i]/CA->Y[i]),2))/2)*CA->n_group[i];
       }
       //aod->addConst -= tmp;
    } else {
      //TODO
      GOFanal->disttype = distribution::log_normal;  //previous GOFanal is always normal distribution
      bmdsConvertSStat(CA, GOFanal, false);  //recalculate using with lognormal transformation
      double divisorTerm = GOFanal->Y[0];
      //rescale to match BMDS 3.x
      for (int i=0; i<GOFanal->n; i++){
         GOFanal->Y[i] /= divisorTerm;
         GOFanal->sd[i] /= divisorTerm;
      }
      //convert from logscale
      double tmpMean, tmpSD;
      for (int i=0; i<GOFanal->n; i++){
        tmpMean = GOFanal->Y[i];
        tmpSD = GOFanal->sd[i];
        GOFanal->Y[i] = exp(tmpMean+pow(tmpSD,2)/2);
        GOFanal->sd[i] = GOFanal->Y[i]*sqrt(exp(pow(tmpSD,2))-1);  

        //reverse response scaling
        GOFanal->Y[i] *= divisorTerm;
        GOFanal->sd[i] *= divisorTerm;

        //convert to log scale
        tmpMean = GOFanal->Y[i];
        tmpSD = GOFanal->sd[i];
        GOFanal->Y[i] = log(tmpMean) - log(1+pow(tmpSD/tmpMean,2))/2;
        GOFanal->sd[i] = sqrt(log(1+pow(tmpSD/tmpMean,2))); 

        tmp += GOFanal->Y[i] * GOFanal->n_group[i];
      }

      //aod->addConst -= tmp;

    }
    aod->addConst -= tmp;
  }

  calcTestsOfInterest(aod);
}



void calcTestsOfInterest(struct continuous_AOD *aod){

  double dev;
  int df;

  //Test #1 - A2 vs Reduced - does mean and/or variance differ across dose groups
  //TOI->llRatio[0] = dev = 2 * (aod->LL[1] - aod->LL[4]);
  dev = (aod->LL[1] - aod->LL[4]);
  //handle underflow/negative zero
  if (dev < BMDS_EPS){
     dev = 0.0;
  }
  aod->TOI.llRatio[0] = dev = 2*dev;

  aod->TOI.DF[0] = df = aod->nParms[1] - aod->nParms[4];
  aod->TOI.pVal[0] = (std::isnan(dev) || dev < 0.0 || df < 0.0) ? -1 : 1.0 - gsl_cdf_chisq_P(dev, df);

  //Test #2 - A1 vs A2 - homogeneity of variance across dose groups
  dev = (aod->LL[1] - aod->LL[0]);
  //handle underflow/negative zero
  if (dev < BMDS_EPS){
     dev = 0.0;
  }
  aod->TOI.llRatio[1] = dev = 2*dev;
  //TOI->llRatio[1] = dev = 2 * (aod->LL[1] - aod->LL[0]);
  aod->TOI.DF[1] = df = aod->nParms[1] - aod->nParms[0];
  aod->TOI.pVal[1] = (std::isnan(dev) || dev < 0.0 || df < 0.0) ? -1 : 1.0 - gsl_cdf_chisq_P(dev, df);

  //Test #3 - A2 vs A3 - Does the model describe variances adequately
  dev = (aod->LL[1] - aod->LL[2]);
  //handle underflow/negative zero
  if (dev < BMDS_EPS){
     dev = 0.0;
  }
  //TOI->llRatio[2] = dev = 2 * (aod->LL[1] - aod->LL[2]);
  aod->TOI.llRatio[2] = dev = 2*dev;
  aod->TOI.DF[2] = df = aod->nParms[1] - aod->nParms[2];
  aod->TOI.pVal[2] = (std::isnan(dev) || dev < 0.0 || df < 0.0) ? -1 : 1.0 - gsl_cdf_chisq_P(dev, df);

  //Test #4 - A3 vs Fitted - does the fitted model describe the obs data adequately 
  dev = (aod->LL[2] - aod->LL[3]);
  //handle underflow/negative zero
  if (dev < BMDS_EPS){
     dev = 0.0;
  }
  //TOI->llRatio[3] = dev = 2 * (aod->LL[2] - aod->LL[3]);
  aod->TOI.llRatio[3] = dev = 2*dev;
  aod->TOI.DF[3] = df = aod->nParms[2] - aod->nParms[3];
  aod->TOI.pVal[3] = (std::isnan(dev) || dev < 0.0 || df < 0.0) ? -1 : 1.0 - gsl_cdf_chisq_P(dev, df);

  //DOF check for test 4
  if (aod->TOI.DF[3] <= 0) {
     aod->TOI.pVal[3] = BMDS_MISSING;
  }
}

void determineAdvDir(struct continuous_analysis *CA){
  int n_rows;

  struct continuous_analysis CAnew;
  //allocate memory for individual size.  Will not use entire array for summary data.
  double* tmpD = (double*)malloc(CA->n * sizeof(double));
  double* tmpY = (double*)malloc(CA->n * sizeof(double));
  double* tmpN = (double*)malloc(CA->n * sizeof(double));
  double* tmpSD = (double*)malloc(CA->n * sizeof(double));
 
  CAnew.doses = tmpD;
  CAnew.Y = tmpY;
  CAnew.n_group = tmpN;
  CAnew.sd = tmpSD;

  if(!CA->suff_stat){
    bmdsConvertSStat(CA, &CAnew,true);
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
  
  free(tmpD);
  free(tmpY);
  free(tmpN);
  free(tmpSD);
}



void bmdsConvertSStat(struct continuous_analysis *CA, struct continuous_analysis *CAss, bool clean){
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
    bool isLognormal = CAss->disttype == distribution::log_normal;
    if (canConvert){

      //bool isLognormal = CA->disttype == distribution::log_normal;
      if (clean){
        if (isLognormal) {
          //Y = cleanSuffStat(SSTAT,X,false);
          Y = cleanSuffStat(SSTAT_LN,X,true,false);
        } else {
//        //Y = cleanSuffStat(SSTAT_LN,X,true);
          Y = cleanSuffStat(SSTAT,X,false,false);
        }
      } else {
        if (isLognormal) {
           Y = SSTAT_LN;
        } else {
           Y = SSTAT;
        }
      }
      n_rows = X.rows();

      for(int i=0; i<n_rows; i++){
        CAss->doses[i] = X(i);
        CAss->Y[i] = Y.col(0)[i];
        CAss->n_group[i] = Y.col(1)[i];
        CAss->sd[i] = Y.col(2)[i];
      }
      CAss->n = n_rows;
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


//replace any double values containing NaN or infinity with BMDS_MISSING
void clean_dicho_results(struct dichotomous_model_result *res, struct dichotomous_GOF *gof, struct BMDS_results *bmdsRes, struct dicho_AOD *aod) {

  //dichotomous_model_result
  for (int i=0; i<res->nparms; i++){
    cleanDouble(&res->parms[i]);
  }   

  for (int i=0; i<res->nparms*res->nparms; i++){
    cleanDouble(&res->cov[i]);
  }

  cleanDouble(&res->max);
  cleanDouble(&res->model_df);
  cleanDouble(&res->total_df);
  cleanDouble(&res->bmd);

  for (int i=0; i<res->dist_numE*2; i++){
    cleanDouble(&res->bmd_dist[i]);
  }

  //dichotomous_GOF
  for (int i=0; i<gof->n; i++){
    cleanDouble(&gof->expected[i]);
//    cleanDouble(&gof->residual[i]);
    cleanDouble(&gof->residual[i]);
    cleanDouble(&gof->ebLower[i]);
    cleanDouble(&gof->ebUpper[i]);
  }
  cleanDouble(&gof->test_statistic);
  cleanDouble(&gof->p_value);
  cleanDouble(&gof->df);
 
  //BMDS_results
  cleanDouble(&bmdsRes->BMD);
  cleanDouble(&bmdsRes->BMDL);
  cleanDouble(&bmdsRes->BMDU);
  cleanDouble(&bmdsRes->AIC);
  cleanDouble(&bmdsRes->BIC_equiv);
  cleanDouble(&bmdsRes->chisq);

  for (int i=0; i<res->nparms; i++){
    cleanDouble(&bmdsRes->stdErr[i]);
    cleanDouble(&bmdsRes->lowerConf[i]);
    cleanDouble(&bmdsRes->upperConf[i]);
  }

  //dicho_AOD
  cleanDouble(&aod->fullLL);
  cleanDouble(&aod->redLL);
  cleanDouble(&aod->fittedLL);
  cleanDouble(&aod->devFit);
  cleanDouble(&aod->devRed);
  cleanDouble(&aod->pvFit);
  cleanDouble(&aod->pvRed);

}

void clean_cont_results(struct continuous_model_result *res, struct BMDS_results *bmdsRes, struct continuous_AOD *aod, struct continuous_GOF *gof){

  //continuous_model_result
  for (int i=0; i<res->nparms; i++){
    cleanDouble(&res->parms[i]);
  }
  for (int i=0; i<res->nparms*res->nparms; i++){
    cleanDouble(&res->cov[i]);
  } 

  cleanDouble(&res->max);
  cleanDouble(&res->model_df);
  cleanDouble(&res->total_df);
  cleanDouble(&res->bmd);
  for (int i=0; i<res->dist_numE*2; i++){
    cleanDouble(&res->bmd_dist[i]);
  }

  //BMDS_results
  cleanDouble(&bmdsRes->BMD);
  cleanDouble(&bmdsRes->BMDL);
  cleanDouble(&bmdsRes->BMDU);
  cleanDouble(&bmdsRes->AIC);
  cleanDouble(&bmdsRes->BIC_equiv);
  cleanDouble(&bmdsRes->chisq);
  for (int i=0; i<res->nparms; i++){
    cleanDouble(&bmdsRes->stdErr[i]);
    cleanDouble(&bmdsRes->lowerConf[i]);
    cleanDouble(&bmdsRes->upperConf[i]);
  } 

  //continuous_AOD and testsOfInterest
  for (int i=0; i<5; i++){
    cleanDouble(&aod->LL[i]);
    cleanDouble(&aod->AIC[i]);
    for (int j=0; j<4; j++){
      cleanDouble(&aod->TOI.llRatio[j]);
      cleanDouble(&aod->TOI.DF[j]);
      cleanDouble(&aod->TOI.pVal[j]);
    }        
  }
  cleanDouble(&aod->addConst);

  //continuous_GOF
  for (int i=0; i<gof->n; i++){
    cleanDouble(&gof->dose[i]);
    cleanDouble(&gof->size[i]);
    cleanDouble(&gof->estMean[i]);
    cleanDouble(&gof->calcMean[i]);
    cleanDouble(&gof->obsMean[i]);
    cleanDouble(&gof->estSD[i]);
    cleanDouble(&gof->calcSD[i]);
    cleanDouble(&gof->obsSD[i]);
    cleanDouble(&gof->res[i]);
    cleanDouble(&gof->ebLower[i]);
    cleanDouble(&gof->ebUpper[i]);
  }

}

void clean_dicho_MA_results(struct dichotomousMA_result *res, struct BMDSMA_results *bmdsRes) {

  //dichotomousMA_result
  for (int i=0; i<res->nmodels; i++){
    cleanDouble(&res->post_probs[i]);
    for (int j=0; j<res->models[i]->nparms; j++){
      cleanDouble(&res->models[i]->parms[j]);
    }
    for (int j=0; j<res->models[i]->nparms*res->models[i]->nparms; j++){
      cleanDouble(&res->models[i]->cov[j]);
    }
    cleanDouble(&res->models[i]->max);
    cleanDouble(&res->models[i]->model_df);
    cleanDouble(&res->models[i]->total_df);
    cleanDouble(&res->models[i]->bmd);
    for (int j=0; j<res->models[i]->dist_numE*2; j++){
      cleanDouble(&res->models[i]->bmd_dist[j]);
    }
  }
  for (int i=0; i<res->dist_numE*2; i++){
    cleanDouble(&res->bmd_dist[i]);
  }

  //BMDSMA_results
  cleanDouble(&bmdsRes->BMD_MA);
  cleanDouble(&bmdsRes->BMDL_MA);
  cleanDouble(&bmdsRes->BMDU_MA);
  for (int i=0; i<res->nmodels; i++){
    cleanDouble(&bmdsRes->BMD[i]);
    cleanDouble(&bmdsRes->BMDL[i]);
    cleanDouble(&bmdsRes->BMDU[i]);
  }

}



//replace any double values containing NaN or infinity with BMDS_MISSING
void cleanDouble(double *val){
  if(!isfinite(*val)) {
     *val = BMDS_MISSING;
  }
}


string BMDS_ENTRY_API __stdcall version(){
  return BMDS_VERSION;   
}


int BMDS_ENTRY_API __stdcall add2(int i, int j) {
    return i + j;
}

void BMDS_ENTRY_API __stdcall testFun(struct test_struct *t){
    std::cout<<"Inside testFun"<<std::endl;

    std::cout<<"t->BMD="<<t->BMD<<std::endl;
    std::cout<<"t->n="<<t->n<<std::endl;
    std::cout<<"t->validResult="<<(t->validResult ? "true" : "false")<<std::endl;
    std::cout<<"Doses:"<<std::endl;
    for (double dose: t->doses) std::cout<< dose << std::endl;
    std::cout<<"Finish testFun"<<std::endl;

}


void convertToPythonDichoRes(struct dichotomous_model_result *res, struct python_dichotomous_model_result *pyRes){
  
  pyRes->model = res->model;
  pyRes->nparms = res->nparms;
  pyRes->parms.assign(res->parms, res->parms + res->nparms);
  pyRes->cov.assign(res->cov, res->cov + res->nparms*res->nparms);
  pyRes->max = res->max;
  pyRes->dist_numE = res->dist_numE;
  pyRes->model_df = res->model_df;
  pyRes->total_df = res->total_df;
  pyRes->bmd_dist.assign(res->bmd_dist, res->bmd_dist + res->dist_numE*2);
  pyRes->bmd = res->bmd;
  pyRes->gof_p_value = res->gof_p_value;
  pyRes->gof_chi_sqr_statistic = res->gof_chi_sqr_statistic;
  
}

void convertFromPythonDichoAnalysis(struct dichotomous_analysis *anal, struct python_dichotomous_analysis *pyAnal){

  anal->model = pyAnal->model;
  anal->n = pyAnal->n;
  anal->BMD_type = pyAnal->BMD_type;
  anal->BMR = pyAnal->BMR;
  anal->alpha = pyAnal->alpha;
  anal->degree = pyAnal->degree;
  anal->samples = pyAnal->samples;
  anal->burnin = pyAnal->burnin;
  anal->parms = pyAnal->parms;
  anal->prior_cols = pyAnal->prior_cols;

  if(pyAnal->n == pyAnal->doses.size() && pyAnal->doses.size() == pyAnal->Y.size() && pyAnal->doses.size() == pyAnal->n_group.size()){
    for (int i=0; i<pyAnal->n; i++){
      anal->Y[i] = pyAnal->Y[i];
      anal->doses[i] = pyAnal->doses[i];
      anal->n_group[i] = pyAnal->n_group[i];
    }
  }
  if(pyAnal->prior.size() > 0){
    for (int i=0; i<pyAnal->prior.size(); i++){
      anal->prior[i] = pyAnal->prior[i];
    }
  }
}

void convertFromPythonDichoRes(struct dichotomous_model_result *res, struct python_dichotomous_model_result *pyRes){
  res->model = pyRes->model;
  res->nparms = pyRes->nparms;
  res->max = pyRes->max;
  res->dist_numE = pyRes->dist_numE;
  res->model_df = pyRes->model_df;
  res->total_df = pyRes->total_df;
  res->bmd = pyRes->bmd;
  res->gof_p_value = pyRes->gof_p_value;
  res->gof_chi_sqr_statistic = pyRes->gof_chi_sqr_statistic;
  //arrays & vectors
  if(pyRes->parms.size() > 0){
    for (int i=0; i<pyRes->parms.size(); i++){
      res->parms[i] = pyRes->parms[i];
    } 
  }
  if(pyRes->cov.size() > 0){
    for (int i=0; i<pyRes->cov.size(); i++){
      res->cov[i] =  pyRes->cov[i];
    }
  }
  if(pyRes->bmd_dist.size() > 0){
    for (int i=0; i<pyRes->bmd_dist.size(); i++){
      res->bmd_dist[i] = pyRes->bmd_dist[i];
    }
  }
}


//void convertFromPythonNestedAnalysis(struct nested_analysis *anal, struct python_nested_analysis *pyAnal){
//
//  int n = pyAnal->doses.size();
//  anal->model = pyAnal->model;
//  anal->n = n;
//  anal->BMD_type = pyAnal->BMD_type;
//  anal->BMR = pyAnal->BMR;
//  anal->LSC_type = pyAnal->LSC_type;
//  anal->ILC_type = pyAnal->ILC_type;
//  anal->alpha = pyAnal->alpha;
//  anal->parms = pyAnal->parms;
//  anal->prior_cols = pyAnal->prior_cols;
//
//  if(n == pyAnal->incidence.size() && n == pyAnal->litterSize.size() && n == pyAnal->lsc.size()){
//    for (int i=0; i<n; i++){
//      anal->Y[i] = pyAnal->incidence[i];
//      anal->doses[i] = pyAnal->doses[i];
//      anal->n_litter[i] = pyAnal->litterSize[i];
//      anal->lsc[i] = pyAnal->lsc[i];
//    }
//  }
//  if(pyAnal->prior.size() > 0){
//    for (int i=0; i<pyAnal->prior.size(); i++){
//      anal->prior[i] = pyAnal->prior[i];
//    }
//  }
//}
//
//
//void convertFromPythonNestedRes(struct nested_model_result *res, struct python_nested_result *pyRes){
//  res->model = pyRes->model;
//  res->nparms = pyRes->nparms;
//  res->max = pyRes->max;
//  res->dist_numE = pyRes->dist_numE;
//  res->model_df = pyRes->model_df;
//  res->total_df = pyRes->total_df;
//  res->bmd = pyRes->bmd;
//  res->comb_p_value = pyRes->combPVal;
//  res->obs_chi_sqr = pyRes->obsChiSq;
//  res->fixedLSC = pyRes->fixedLSC;
//  //arrays & vectors
//  if(pyRes->parms.size() > 0){
//    for (int i=0; i<pyRes->parms.size(); i++){
//      res->parms[i] = pyRes->parms[i];
//    } 
//  }
//  if(pyRes->cov.size() > 0){
//    for (int i=0; i<pyRes->cov.size(); i++){
//      res->cov[i] =  pyRes->cov[i];
//    }
//  }
//  if(pyRes->bmd_dist.size() > 0){
//    for (int i=0; i<pyRes->bmd_dist.size(); i++){
//      res->bmd_dist[i] = pyRes->bmd_dist[i];
//    }
//  }
//}

void convertFromPythonDichoMAAnalysis(struct dichotomousMA_analysis *MA, struct python_dichotomousMA_analysis *pyMA){

  MA->nmodels = pyMA->nmodels;
  for (int i=0; i<pyMA->nmodels; i++){
    MA->priors[i] = pyMA->priors[i].data();
  } 
  MA->nparms = pyMA->nparms.data();
  MA->actual_parms = pyMA->actual_parms.data();
  MA->prior_cols = pyMA->prior_cols.data();
  MA->models = pyMA->models.data();
  MA->modelPriors = pyMA->modelPriors.data();
 

}

void convertFromPythonDichoMARes(struct dichotomousMA_result *res, struct python_dichotomousMA_result *pyRes){

  res->nmodels = pyRes->nmodels;
  for (int i=0; i<pyRes->nmodels; i++){
    //res->models[i] = pyRes->models[i];
    convertFromPythonDichoRes(res->models[i], &pyRes->models[i]);
  }
  res->dist_numE = pyRes->dist_numE;
  res->post_probs = pyRes->post_probs.data();
  res->bmd_dist = pyRes->bmd_dist.data();

}

void convertFromPythonContAnalysis(struct continuous_analysis *anal, struct python_continuous_analysis *pyAnal){

  anal->model = pyAnal->model;
  anal->n = pyAnal->n;
  anal->BMD_type = pyAnal->BMD_type;
  anal->isIncreasing;
  anal->BMR = pyAnal->BMR;
  anal->tail_prob = pyAnal->tail_prob;
  anal->disttype = pyAnal->disttype;
  anal->alpha = pyAnal->alpha;
  anal->degree = pyAnal->degree;
  anal->samples = pyAnal->samples;
  anal->burnin = pyAnal->burnin;
  anal->parms = pyAnal->parms;
  anal->prior_cols = pyAnal->prior_cols;
  anal->transform_dose = pyAnal->transform_dose;
  anal->suff_stat = pyAnal->suff_stat;

  if(pyAnal->n == pyAnal->doses.size() && pyAnal->doses.size() == pyAnal->Y.size() && pyAnal->doses.size() == pyAnal->n_group.size()){
    for (int i=0; i<pyAnal->n; i++){
      anal->Y[i] = pyAnal->Y[i];
      anal->doses[i] = pyAnal->doses[i];
      anal->n_group[i] = pyAnal->n_group[i];
      anal->sd[i] = pyAnal->sd[i];
    }
  }
  if(pyAnal->prior.size() > 0){
    for (int i=0; i<pyAnal->prior.size(); i++){
      anal->prior[i] = pyAnal->prior[i];
    }
  }
}

void convertFromPythonContRes(struct continuous_model_result *res, struct python_continuous_model_result *pyRes){
  res->model = pyRes->model;
  res->dist = pyRes->dist;
  res->nparms = pyRes->nparms;
  res->max = pyRes->max;
  res->dist_numE = pyRes->dist_numE;
  res->model_df = pyRes->model_df;
  res->total_df = pyRes->total_df;
  res->bmd = pyRes->bmd;
  //arrays & vectors
  if(pyRes->parms.size() > 0){
    for (int i=0; i<pyRes->parms.size(); i++){
      res->parms[i] = pyRes->parms[i];
    } 
  }
  if(pyRes->cov.size() > 0){
    for (int i=0; i<pyRes->cov.size(); i++){
      res->cov[i] =  pyRes->cov[i];
    }
  }
  if(pyRes->bmd_dist.size() > 0){
    for (int i=0; i<pyRes->bmd_dist.size(); i++){
      res->bmd_dist[i] = pyRes->bmd_dist[i];
    }
  }
}

void convertToPythonContRes(struct continuous_model_result *res, struct python_continuous_model_result *pyRes){
  
  pyRes->model = res->model;
  pyRes->dist = res->dist;
  pyRes->nparms = res->nparms;
  pyRes->parms.assign(res->parms, res->parms + res->nparms);
  pyRes->cov.assign(res->cov, res->cov + res->nparms*res->nparms);
  pyRes->max = res->max;
  pyRes->dist_numE = res->dist_numE;
  pyRes->model_df = res->model_df;
  pyRes->total_df = res->total_df;
  pyRes->bmd_dist.assign(res->bmd_dist, res->bmd_dist + res->dist_numE*2);
  pyRes->bmd = res->bmd;
  
}

void BMDS_ENTRY_API __stdcall pythonBMDSDicho(struct python_dichotomous_analysis *pyAnal, struct python_dichotomous_model_result *pyRes){

  //1st convert from python struct
  dichotomous_analysis anal;
  anal.Y = new double[pyAnal->n];
  anal.doses = new double[pyAnal->n];
  anal.n_group = new double[pyAnal->n];
  anal.prior = new double[pyAnal->parms*pyAnal->prior_cols];
  convertFromPythonDichoAnalysis(&anal, pyAnal); 

  dichotomous_model_result res;
  res.parms = new double[pyRes->nparms];
  res.cov = new double[pyRes->nparms*pyRes->nparms];
  res.bmd_dist = new double[pyRes->dist_numE*2];
  convertFromPythonDichoRes(&res, pyRes);

  runBMDSDichoAnalysis(&anal, &res, &pyRes->gof, &pyRes->bmdsRes, &pyRes->aod);     

  convertToPythonDichoRes(&res, pyRes);

}

void BMDS_ENTRY_API __stdcall pythonBMDSDichoMA(struct python_dichotomousMA_analysis *pyMA, struct python_dichotomousMA_result *pyRes){

//convert python_dichtomousMA_analysis to dichotomousMA_analysis
  dichotomousMA_analysis MA;
  MA.priors = new double *[pyMA->nmodels];
  MA.nparms = new int[pyMA->nmodels];
  MA.actual_parms = new int[pyMA->nmodels];
  MA.prior_cols = new int[pyMA->nmodels];
  MA.models = new int[pyMA->nmodels];
  MA.modelPriors = new double[pyMA->nmodels];
  convertFromPythonDichoMAAnalysis(&MA, pyMA);

//convert python_dichotomous_analysis to dichotomous_analysis
  dichotomous_analysis DA;
  DA.Y = new double[pyMA->pyDA.n];
  DA.doses = new double[pyMA->pyDA.n];
  DA.n_group = new double[pyMA->pyDA.n];
  DA.prior = new double[pyMA->pyDA.parms*pyMA->pyDA.prior_cols];
  convertFromPythonDichoAnalysis(&DA, &pyMA->pyDA);

//  convertFromPythonDichoMARes(&res, pyRes);

  struct dichotomous_model_result** res = new struct dichotomous_model_result* [pyRes->nmodels];
  for (int i=0; i<pyRes->nmodels; i++){
    res[i] = (struct dichotomous_model_result*)malloc(sizeof(struct dichotomous_model_result));
    res[i]->model = pyRes->models[i].model;
    res[i]->nparms = pyRes->models[i].nparms;
    res[i]->dist_numE = pyRes->models[i].dist_numE;
    res[i]->parms = (double*)malloc(sizeof(double)*pyRes->models[i].nparms);
    res[i]->cov = (double*)malloc(sizeof(double)*pyRes->models[i].nparms*pyRes->models[i].nparms);
    res[i]->bmd_dist = (double*)malloc(sizeof(double)*pyRes->models[i].dist_numE*2);
  }

  struct dichotomousMA_result maRes;
  maRes.nmodels = pyRes->nmodels;
  maRes.models = res;
  maRes.dist_numE = pyRes->dist_numE;
  maRes.post_probs = new double[pyRes->nmodels];
  maRes.bmd_dist = new double[pyRes->dist_numE*2];

  runBMDSDichoMA(&MA, &DA, &maRes, &pyRes->bmdsRes);
//convert back to python objs
//pyDA should not change
  pyRes->post_probs.resize(pyRes->nmodels);
  for(int i=0; i<pyRes->nmodels; i++){
    pyRes->post_probs[i] = maRes.post_probs[i];
    pyRes->models[i].max = maRes.models[i]->max;
    pyRes->models[i].model_df = maRes.models[i]->model_df;
    pyRes->models[i].total_df = maRes.models[i]->total_df;
    pyRes->models[i].bmd = maRes.models[i]->bmd;                  
    pyRes->models[i].gof_p_value = maRes.models[i]->gof_p_value;           
    pyRes->models[i].gof_chi_sqr_statistic = maRes.models[i]->gof_chi_sqr_statistic;  
    int nparms = pyRes->models[i].nparms;

    pyRes->models[i].parms.resize(nparms);
    for(int j=0; j<nparms; j++){
      pyRes->models[i].parms[j] = maRes.models[i]->parms[j];
    }
    pyRes->models[i].cov.resize(nparms*nparms);
    for(int j=0; j<nparms*nparms; j++){
      pyRes->models[i].cov[j] = maRes.models[i]->cov[j];
    }
    pyRes->models[i].bmd_dist.resize(pyRes->dist_numE*2);
    for(int j=0; j<pyRes->dist_numE*2; j++){
      pyRes->models[i].bmd_dist[j] = maRes.models[i]->bmd_dist[j];
    } 
  } 
  pyRes->bmd_dist.resize(pyRes->dist_numE*2);
  for(int i=0; i<pyRes->dist_numE*2; i++){
    pyRes->bmd_dist[i] = maRes.bmd_dist[i];
  }

}

void BMDS_ENTRY_API __stdcall pythonBMDSCont(struct python_continuous_analysis *pyAnal, struct python_continuous_model_result *pyRes){

  //convert pyAnal to anal
  continuous_analysis anal;
  anal.Y = new double[pyAnal->n];
  anal.doses = new double[pyAnal->n];
  anal.sd = new double[pyAnal->n];
  anal.n_group = new double[pyAnal->n];
  anal.prior = new double[pyAnal->parms*pyAnal->prior_cols];
  convertFromPythonContAnalysis(&anal, pyAnal);

  //convert pyRes to res
  continuous_model_result res;
  res.parms = new double[pyRes->nparms];
  res.cov = new double[pyRes->nparms*pyRes->nparms];
  res.bmd_dist = new double[pyRes->dist_numE*2];
  convertFromPythonContRes(&res, pyRes);

  runBMDSContAnalysis(&anal, &res, &pyRes->bmdsRes, &pyRes->aod, &pyRes->gof, &pyAnal->detectAdvDir, &pyAnal->restricted);

  convertToPythonContRes(&res, pyRes);

}

void BMDS_ENTRY_API __stdcall pythonBMDSMultitumor(struct python_multitumor_analysis *pyAnal, struct python_multitumor_result *pyRes){
 

   //run each individual multistage model
   for (int i=0;i<pyAnal->ndatasets;i++){
     for (int j=0; j<pyAnal->nmodels[i]; j++){
       pythonBMDSDicho(&pyAnal->models[i][j], &pyRes->models[i][j]);  
     } 
   }
   
   //select models
   selectMultitumorModel();
   //dev result stubs
   pyRes->selectedModelIndex.push_back(0);
   pyRes->selectedModelIndex.push_back(2);
   pyRes->selectedModelIndex.push_back(2);

   //run MSCombo
   runMultitumorModel();  
   //dev result stubs
   pyRes->BMD = 6.005009;
   pyRes->BMDL = 5.402722;
   pyRes->BMDU = 6.9757479;
   pyRes->slopeFactor = 0.018509189;
   pyRes->combined_LL = -563.8220487; 
   pyRes->combined_LL_const = 521.8050805; 
  
}  


void selectMultitumorModel(){
   std::cout<<"selectMultitumorModel to be implemented"<<std::endl;
}

void runMultitumorModel(){
   std::cout<<"runMultitumorModel"<<std::endl;
}


double PROBABILITY_INRANGE(double ex){

  if(ex <= 0.00000001){
    ex = 1.0e-7;
  }
  if(ex >= 0.9999999){
    ex = 0.9999999;
  }
  return ex;
}

void BMDS_ENTRY_API __stdcall pythonBMDSNested(struct python_nested_analysis *pyAnal, struct python_nested_result *pyRes){



   //////////////////////////////
   // Parm vector indexes
   // 0 = alpha (intercept or background)
   // 1 = beta
   // 2 = THETA1
   // 3 = THETA2
   // 4 = rho
   // 5 = PHISTART - first PHI parameter
   //
   int BGINDEX = 0;  //location of background parameter in parms vector
   int PHISTART = 6; //location of first PHI parameter in parms vector
   int THETA1_INDEX = 3; //location of THETA1 parameter in parms vector 
   int THETA2_INDEX = 4; //location of THETA2 parameter in parms vector
   bool validResult = false;
 
   //set seed from time clock if default seed=0 is specified
   if (pyAnal->seed == BMDS_MISSING){
     pyAnal->seed = time (NULL);
   }


   int Nobs = pyAnal->doses.size();
   std::vector<double> Xi(Nobs); //doses (independent data var)
   std::vector<double> Yp(Nobs); //incidence (positive dependent var)
   std::vector<double> Yn(Nobs); //negative dependent var
   std::vector<double> Ls(Nobs); //litter size
   std::vector<double> SR(Nobs); //scaled residual
   std::vector<int> Xg(Nobs, BMDS_MISSING); //markings for group number

   //Need to implement
   //quick workaround assuming vectors are already sorted
   Xi = pyAnal->doses;
   Yp = pyAnal->incidence;
   for (int i=0; i<Nobs; i++){
     Yn[i] = pyAnal->litterSize[i] -  pyAnal->incidence[i];
   }
   Ls = pyAnal->litterSize;
   //Sort_4_By_Dose(Nobs, Xi, Yn, Yp, Ls);


   //get correct number of dose groups
   //assumes vectors are already sorted
//   int ngrp = 0;  //number of litters in dose group
//   double tmp = BMDS_MISSING;
//   int groupnum = 0;  //group number
//   std::vector<int> grpSize;
//   for (int i=0; i<Nobs; i++){
//      if (Xi[i] == tmp){
//        //existing dose group
//        ngrp++;
//      } else {
//        //new dose group
//        if (i!=0){
//          grpSize.push_back(ngrp);
//        }
//        tmp = Xi[i];
//        ngrp = 1;
//        groupnum++;
//      }
//      Xg[i] = groupnum;
//   }
   double junk1 = Xi[0];
   int junk = 1;
   for (int i=0; i<Nobs; i++){
      if (Xi[i] == junk1){
	 Xg[i] = junk;
      } else {
	 Xg[i] = ++junk;
	 junk1 = Xi[i];
      }
   }
   int ngrp=junk;

   std::cout<<"ngrp = "<<ngrp<<std::endl;
   //push last group size;
   //grpSize.push_back(ngrp);
   //pyRes->nparms = 5 + grpSize.size();
   pyRes->nparms = 5 + ngrp;

   std::cout<<"Xi   , Yp,   Yn,   Ls,   Xg"<<std::endl;
   for (int i=0; i<Nobs; i++){
      std::cout<<Xi[i]<<", "<<Yp[i]<<", "<<Yn[i]<<", "<<Ls[i]<<", "<<Xg[i]<<std::endl;
   }

//   int Nobs = 0;
//   for (int i=0; i<grpSize.size(); i++){
//      Nobs += grpSize[i];
//   }

   //create and initialize vectors needed
   pyRes->parms.resize(pyRes->nparms);
   for (int i=0; i<pyRes->nparms; i++){
     pyRes->parms[i] = BMDS_MISSING;
   }
//   std::vector<bool> Spec(pyRes->nparms, false);
//   std::vector<bool> iniP(pyRes->nparms, false);
//   std::vector<VarList> varsum(3);
//   std::vector<AnaList> anasum(3);
//   std::vector<double> rlevel(5);
//   std::vector<double> bmdl(5);
//   std::vector<std::vector<double>> vcv(pyRes->nparms, std::vector<double>(pyRes->nparms));
//   double sijfixed;
//   std::vector<bool> SpBak;  //old values for spec
//   std::vector<double> pBak;  //old values for parms[];
//
//
//   std::cout<<"numGrps = "<<grpSize.size()<<std::endl;
//   for(int i=0; i<grpSize.size(); i++){
//      std::cout<<"i:"<<i<<", grpSize:"<<grpSize[i]<<std::endl;
//   }
//   std::cout<<"seed="<<pyAnal->seed<<std::endl;
  
  int knownParms = 0;
  //handle specified parms
  if (!pyAnal->estBackground){
    //set background to zero
    pyRes->parms[BGINDEX] = 0;
//    Spec[BGINDEX] = true;
//    iniP[BGINDEX] = true;  //Specified implies initialized
    knownParms++;
  }   
  if (pyAnal->ILC_type != 1){
    //set phi parm values to zero
    for (int i=PHISTART; i<pyRes->nparms; i++){
      pyRes->parms[i] = 0.0;
//      Spec[i] = true;
//      iniP[i] = true;
    }
    knownParms += ngrp;
  }
  if (!pyAnal->useLSC){
    //set theta parm values to zero
    for (int i=THETA1_INDEX; i<=THETA2_INDEX; i++){
       pyRes->parms[i] = 0.0;
//       Spec[i] = true;
//       iniP[i] = true;
    }
    knownParms += 2;  //theta1 and theta2
  }
   
  if(Nobs < pyRes->nparms - knownParms){
    std::cout<<"Error: Fewer observations "<<Nobs<<" than estimated parameters "<<pyRes->nparms-knownParms<<std::endl;
    return;
  }

  //make sure all litter sizes > 0
  for (int i=0; i<Nobs; i++){
    if (pyAnal->litterSize[i] <= 0.0){
      std::cout<<"Error: All litter sizes must be greater than zero."<<std::endl;
      return;
    }
  }

  //compute init_lkf for full model and init_lkr for reduced model
  double lkf = 0.0;
  double pdep = 0.0;  //positive response
  double ndep = 0.0;  //negative response 
//  varsum[0].S = 0;
//  varsum[1].S = 0;
  double W;  
  for (int i=0; i<Nobs; i++){
//    varsum[0].S += Yp[i];
//    varsum[1].S += Yn[i];
    pdep += Yp[i];
    ndep += Yn[i];
    W = Yp[i] / (Yp[i] + Yn[i]);
    if (W > 0){
      lkf += Yp[i] * log (W);
    }
    if (W < 1){
      lkf += Yn[i] * log(1-W);
    }
  }
//  W = varsum[0].S / (varsum[0].S + varsum[1].S);
  W = pdep / (pdep + ndep);
//  double lkr = varsum[0].S * log(W) + varsum[1].S * log(1-W);
  double lkr = pdep * log(W) + ndep * log(1-W);

  //Nlogist_fit(nparm, ngrp, Parms, EPS, &iter, &xlk);
  double EPS = 3.0e-8;
  int iter = 0;
  double xlk = 0.0;
 // int ret = Nlogist_fit(pyRes->nparms, ngrp, pyRes->parms, EPS, &iter, &xlk, Xg, Ls, Xi);
////////////////////////////////////////////////////
// START NLOGIST_FIT
// //////////////////////////////////////////////////
 
  std::cout<<"Inside Nlogist_fit"<<std::endl;
  double sum1 = 0.0;
  double nij = 0.0;
  int index=0;
  while (Xg[index] == 1){
    sum1 += Ls[index];
    nij +=1.0;
    index++;
  }
  double smean1 = sum1/nij;  //the average lsc in group 1
  std::cout<<"average lsc in group 1:"<<smean1<<std::endl;

  sum1 = 0.0;
  nij = 0.0;
  double smax = Ls[0];
  double smin = Ls[0];
  double xmax = 0.0;
  double x = 0.0;
  for (int i=0; i<=Nobs; i++){
    x = Xi[i];
    sum1 += Ls[i];
    nij += 1.0;
    if (Ls[i] > smax) smax = Ls[i];
    if (Ls[i] < smin) smin = Ls[i];
    if (x > xmax) xmax = x;
  }
  double smean = sum1/nij;  //overall average lsc
  std::cout<<"overall average lsc:"<<smean<<std::endl;

  double sijfixed = smean;  //default to overall mean.  Not used if pyAnal->LSCType = 0 (do not use LSC) 
  if (pyAnal->LSC_type == 2){  //control group mean
     sijfixed = smean1;
  }

  std::vector<double> pBak = pyRes->parms;

  //compute default starting values for all unfixed parameters
  //revisit if we need to fix any parameters

  //Get the starting values in stages
  //First, get starting values for alpha, beta, and rho


   std::vector<double> tmpXi(Nobs); //doses (independent data var)
   std::vector<double> tmpYp(Nobs); //incidence (positive dependent var)
   std::vector<double> tmpYn(Nobs); //negative dependent var
   std::vector<double> tmpy(2);
   //std::vector<std::vector<double>> tmpvcv(2,std::vector<double>(2));
   Eigen::Matrix2d tmpvcv(2,2);

  index = 0;
  //convert nested data to dichotomous data
  for (int j = 0; j<ngrp; j++){
     tmpYn[j] = 0;
     tmpYp[j] = 0;
     while (index < Nobs && Xg[index] == (j+1)){  //note - index is 1 behind group number (index=0 for group 1)
        tmpYn[j] += Yn[index];
	tmpYp[j] += Yp[index];
	tmpXi[j] = Xi[index];
	index++;
     }
     std::cout<<"j:"<<j<<", tmpYn:"<<tmpYn[j]<<", tmpYp:"<<tmpYp[j]<<", tmpXi[j]:"<<tmpXi[j]<<std::endl;
  }
  
  //compute initial estimates
  double ymin = 1.0;
  for (int i=0; i<ngrp; i++){
     W = tmpYp[i]/(tmpYp[i]+tmpYn[i]);
     probability_inrange(&W);
     if (W < ymin) ymin = W;
  }

  for (int i=0; i<ngrp; i++){
     x = tmpXi[i];
     std::cout<<"i:"<<i<<", init x="<<x<<std::endl;
     W = log((tmpYp[i] - ymin*(tmpYp[i] + tmpYn[i]) + 0.5)/(tmpYn[i] + 0.5));
     if (x<=CloseToZero){
        x=Log_zero;
	std::cout<<"x logzero for i="<<i<<std::endl;
     } else {
	x = log(x);
	std::cout<<"x logx for i="<<i<<std::endl;
     }
     //tmpvcv[0][0] += 1.0;
     //tmpvcv[0][1] += x;
     //tmpvcv[1][1] += x*x;
     tmpvcv(0,0) += 1.0;
     tmpvcv(0,1) += x;
     tmpvcv(1,1) += x*x;
     tmpy[0] += W;
     tmpy[1] += W*x;
     std::cout<<"i:"<<i<<", x="<<x<<", W="<<W<<std::endl;
     std::cout<<"tmpvcv(0,0)="<<tmpvcv(0,0)<<", tmpvcv(0,1)="<<tmpvcv(0,1)<<", tmpvcv(1,1)="<<tmpvcv(1,1)<<std::endl;
     std::cout<<"tmpy[0]="<<tmpy[0]<<", tmpy[1]="<<tmpy[1]<<std::endl;
  }
  //tmpvcv[1][0] = tmpvcv[0][1];
  tmpvcv(1,0) = tmpvcv(0,1);
  pyRes->parms[0] = ymin + 0.001;  //in case ymin=0

  std::cout<<"det tmpvcv:"<<tmpvcv.determinant()<<std::endl; 
  if (tmpvcv.determinant() > 0){
    Eigen::Matrix2d invtmpvcv = tmpvcv.inverse();
    pyRes->parms[1] = invtmpvcv(0,0) * tmpy[0] + invtmpvcv(0,1) * tmpy[1];
    pyRes->parms[4] = invtmpvcv(1,0) * tmpy[0] + invtmpvcv(1,1) * tmpy[1];
  } else {
    pyRes->parms[1] = -1.0;
    pyRes->parms[4] = 1.001;
  }

  if (pyAnal->restricted && pyRes->parms[4] < 1.000){
    pyRes->parms[4] = 1.0001;
  }

  //setup log-logistic model
  for (int i=5; i<pyRes->nparms; i++){
     pyRes->parms[i] = 0.0;
  }
  pyRes->parms[2] = pyRes->parms[3] = 0.0;

  //simple non-nested log-logistic model
  double gtol = 3e-8;

  std::cout<<std::endl;
  std::cout<<"b4 Nlogist_lk pass=1"<<std::endl;
  for (int i=0; i<pyRes->nparms; i++){
     std::cout<<"i:"<<i<<", parm:"<<pyRes->parms[i]<<std::endl;
  }

  //This first pass has theta1 p[2] and theta2 p[3] fixed to zero
  int pass = 1;
  double retVal = opt_nlogistic(pyRes->parms, Ls, Xi, Xg, Yp, Yn, smax, smin, pyAnal->restricted, pass);
  std::cout<<"retVal:"<<retVal<<std::endl;


  //Now alpha p[0], beta p[1], and rho[4] contain starting estimates
  //Second, get initial values for Phi's
  int count = 0;
  //for (int i=5; i<n; i++){
  //  count += SpBak[i];
  //}
  //Currently do not allow specification of Phi values, so this is always true
  if (count < ngrp){
    // leave theta1 and theta2 = 0.0;
    for (int i=5; i<pyRes->nparms; i++){
      //set back to original parms
      //pyRes->parms[i] = pBak[i];
      //if (SpBak[i] == 0){
      pyRes->parms[i] = 0.01;
      //}
    }
    for (int i=5; i<pyRes->nparms; i++){
      pyRes->parms[i] = pyRes->parms[i]/(1-pyRes->parms[i]);  //Phi --> Psi
    }
  
    std::cout<<std::endl;
    std::cout<<"b4 Nlogist_lk pass=2"<<std::endl;
    for (int i=0; i<pyRes->nparms; i++){
       std::cout<<"i:"<<i<<", parm:"<<pyRes->parms[i]<<std::endl;
    }

    pass = 2;
    retVal = opt_nlogistic(pyRes->parms, Ls, Xi, Xg, Yp, Yn, smax, smin, pyAnal->restricted, pass);

    //Transform parameters to "external" form
    for (int i=5; i<pyRes->nparms; i++){
      pyRes->parms[i] = pyRes->parms[i]/(1+pyRes->parms[i]);  //Psi --> Phi
    }

    //When theta1 ==0, internal and external forms for alpha are the same,
    //so we do not need to back transform p[0]    
  }

  //Finally, get starting values for Thetas
  for (int i=2; i<=3; i++){
    //if not set
    pyRes->parms[i] = 0.0;
  }

  //Fit the model

  //Transform the parameters to the "internal" form
  for (int i=5; i<pyRes->nparms; i++){
    pyRes->parms[i] = pyRes->parms[i]/(1-pyRes->parms[i]);  //Phi --> Psi
  }

  double junk3;
  //don't currently allow spec so this is always true
//  if (Spec[0] == Spec[2]){
  junk1 = pyRes->parms[0];
  junk3 = pyRes->parms[2];
  pyRes->parms[0] = junk1 + smin * junk3;
  pyRes->parms[2] = junk1 + smax * junk3;
//  }

  //ML fit and return log-likelihood value
  std::cout<<std::endl;
  std::cout<<"b4 Nlogist_lk pass=3"<<std::endl;
  for (int i=0; i<pyRes->nparms; i++){
     std::cout<<"i:"<<i<<", parm:"<<pyRes->parms[i]<<std::endl;
  }
  pass = 3;
  retVal = opt_nlogistic(pyRes->parms, Ls, Xi, Xg, Yp, Yn, smax, smin, pyAnal->restricted, pass);





//  double ll = Nlogist_lk(pyRes->parms, Ls, Xi, Xg, Yp, Yn, smax, smin);

//  std::cout<<"after Nlogist_lk"<<std::endl;
//  for (int i=0; i<pyRes->nparms; i++){
//     std::cout<<"i:"<<i<<", parm:"<<pyRes->parms[i]<<std::endl;
//  }

//   //stubbed results
//   pyRes->bmdsRes.validResult = true;
//   pyRes->bmdsRes.BMD = 12.95166613;
//   pyRes->bmdsRes.BMDL = 9.643478831;
//   pyRes->bmdsRes.BMDU = (double)BMDS_MISSING;
//
//   pyRes->bmdsRes.AIC = 546.9572409;
//   pyRes->bmdsRes.chisq = 19.6087053;
//   pyRes->combPVal = 0.994;
//   //pyRes->df = 35;
//   pyRes->model_df = 35;
//   pyRes->total_df = 35;
//
//   pyRes->nparms = 9;
//   pyRes->parms.push_back(0.084733516);
//   pyRes->parms.push_back(-4.109651594);
//   pyRes->parms.push_back(0.004761366);
//   pyRes->parms.push_back(-0.055489253);
//   pyRes->parms.push_back(1.0);
//   pyRes->parms.push_back(0);
//   pyRes->parms.push_back(0);
//   pyRes->parms.push_back(0);
//   pyRes->parms.push_back(0);
//
//   pyAnal->seed = 1687267999;
//   //pyAnal->iterations = 1000;
//   pyRes->LL = -269.478205;
//   pyRes->obsChiSq = 19.6087053;
//
//   pyRes->boot.numRuns = 3;   
//   pyRes->boot.pVal.push_back(0.994);
//   pyRes->boot.pVal.push_back(0.997);
//   pyRes->boot.pVal.push_back(0.991);
//   pyRes->boot.pVal.push_back(0.994);
//   pyRes->boot.perc50.push_back(38.79787451);
//   pyRes->boot.perc50.push_back(38.19554318);
//   pyRes->boot.perc50.push_back(37.75643018);
//   pyRes->boot.perc50.push_back(38.26776644);
//   pyRes->boot.perc90.push_back(51.1848655);
//   pyRes->boot.perc90.push_back(50.4569082);
//   pyRes->boot.perc90.push_back(50.3883043);
//   pyRes->boot.perc90.push_back(50.6118335);
//   pyRes->boot.perc95.push_back(54.69182);
//   pyRes->boot.perc95.push_back(53.99859);
//   pyRes->boot.perc95.push_back(54.18472);
//   pyRes->boot.perc95.push_back(54.55846);
//   pyRes->boot.perc99.push_back(60.479415);
//   pyRes->boot.perc99.push_back(63.639965);
//   pyRes->boot.perc99.push_back(61.778094);
//   pyRes->boot.perc99.push_back(62.421371);
//
//   pyRes->SRs.push_back(-0.31484);
//   pyRes->SRs.push_back(0.314837);
//   pyRes->SRs.push_back(-0.31484);
//   pyRes->SRs.push_back(0.314837);
//   pyRes->SRs.push_back(-0.31484);
//   pyRes->SRs.push_back(0.314837);
//
//   pyRes->litter.dose.push_back(0);
//   pyRes->litter.dose.push_back(0);
//   pyRes->litter.dose.push_back(0);
//   pyRes->litter.dose.push_back(0);
//   pyRes->litter.dose.push_back(0);
//   pyRes->litter.dose.push_back(0);
//   pyRes->litter.dose.push_back(0);
//   pyRes->litter.dose.push_back(0);
//   pyRes->litter.dose.push_back(0);
//   pyRes->litter.dose.push_back(0);
//   pyRes->litter.dose.push_back(25);
//   pyRes->litter.dose.push_back(25);
//   pyRes->litter.dose.push_back(25);
//   pyRes->litter.dose.push_back(25);
//   pyRes->litter.dose.push_back(25);
//   pyRes->litter.dose.push_back(25);
//   pyRes->litter.dose.push_back(25);
//   pyRes->litter.dose.push_back(25);
//   pyRes->litter.dose.push_back(25);
//   pyRes->litter.dose.push_back(25);
//   pyRes->litter.dose.push_back(50);
//   pyRes->litter.dose.push_back(50);
//   pyRes->litter.dose.push_back(50);
//   pyRes->litter.dose.push_back(50);
//   pyRes->litter.dose.push_back(50);
//   pyRes->litter.dose.push_back(50);
//   pyRes->litter.dose.push_back(50);
//   pyRes->litter.dose.push_back(50);
//   pyRes->litter.dose.push_back(50);
//   pyRes->litter.dose.push_back(50);
//   pyRes->litter.dose.push_back(100);
//   pyRes->litter.dose.push_back(100);
//   pyRes->litter.dose.push_back(100);
//   pyRes->litter.dose.push_back(100);
//   pyRes->litter.dose.push_back(100);
//   pyRes->litter.dose.push_back(100);
//   pyRes->litter.dose.push_back(100);
//   pyRes->litter.dose.push_back(100);
//   pyRes->litter.dose.push_back(100);
//
//   pyRes->litter.LSC.push_back(9);
//   pyRes->litter.LSC.push_back(9);
//   pyRes->litter.LSC.push_back(10);
//   pyRes->litter.LSC.push_back(10);
//   pyRes->litter.LSC.push_back(11);
//   pyRes->litter.LSC.push_back(13);
//   pyRes->litter.LSC.push_back(14);
//   pyRes->litter.LSC.push_back(14);
//   pyRes->litter.LSC.push_back(15);
//   pyRes->litter.LSC.push_back(16);
//   pyRes->litter.LSC.push_back(9);
//   pyRes->litter.LSC.push_back(9);
//   pyRes->litter.LSC.push_back(10);
//   pyRes->litter.LSC.push_back(10);
//   pyRes->litter.LSC.push_back(11);
//   pyRes->litter.LSC.push_back(12);
//   pyRes->litter.LSC.push_back(13);
//   pyRes->litter.LSC.push_back(14);
//   pyRes->litter.LSC.push_back(14);
//   pyRes->litter.LSC.push_back(14);
//   pyRes->litter.LSC.push_back(7);
//   pyRes->litter.LSC.push_back(10);
//   pyRes->litter.LSC.push_back(10);
//   pyRes->litter.LSC.push_back(11);
//   pyRes->litter.LSC.push_back(11);
//   pyRes->litter.LSC.push_back(11);
//   pyRes->litter.LSC.push_back(11);
//   pyRes->litter.LSC.push_back(14);
//   pyRes->litter.LSC.push_back(14);
//   pyRes->litter.LSC.push_back(15);
//   pyRes->litter.LSC.push_back(8);
//   pyRes->litter.LSC.push_back(10);
//   pyRes->litter.LSC.push_back(11);
//   pyRes->litter.LSC.push_back(11);
//   pyRes->litter.LSC.push_back(12);
//   pyRes->litter.LSC.push_back(12);
//   pyRes->litter.LSC.push_back(13);
//   pyRes->litter.LSC.push_back(14);
//   pyRes->litter.LSC.push_back(14);
//
//   pyRes->litter.estProb.push_back(0.127585814);
//   pyRes->litter.estProb.push_back(0.127585814);
//   pyRes->litter.estProb.push_back(0.13234718);
//   pyRes->litter.estProb.push_back(0.13234718);
//   pyRes->litter.estProb.push_back(0.137108547);
//   pyRes->litter.estProb.push_back(0.14663128);
//   pyRes->litter.estProb.push_back(0.151392646);
//   pyRes->litter.estProb.push_back(0.151392646);
//   pyRes->litter.estProb.push_back(0.156154013);
//   pyRes->litter.estProb.push_back(0.160915379);
//   pyRes->litter.estProb.push_back(0.301527034);
//   pyRes->litter.estProb.push_back(0.301527034);
//   pyRes->litter.estProb.push_back(0.297781775);
//   pyRes->litter.estProb.push_back(0.297781775);
//   pyRes->litter.estProb.push_back(0.294373053);
//   pyRes->litter.estProb.push_back(0.291294677);
//   pyRes->litter.estProb.push_back(0.288539894);
//   pyRes->litter.estProb.push_back(0.286101463);
//   pyRes->litter.estProb.push_back(0.286101463);
//   pyRes->litter.estProb.push_back(0.286101463);
//   pyRes->litter.estProb.push_back(0.433391608);
//   pyRes->litter.estProb.push_back(0.410232266);
//   pyRes->litter.estProb.push_back(0.410232266);
//   pyRes->litter.estProb.push_back(0.40315061);
//   pyRes->litter.estProb.push_back(0.40315061);
//   pyRes->litter.estProb.push_back(0.40315061);
//   pyRes->litter.estProb.push_back(0.40315061);
//   pyRes->litter.estProb.push_back(0.38390157);
//   pyRes->litter.estProb.push_back(0.38390157);
//   pyRes->litter.estProb.push_back(0.378161814);
//   pyRes->litter.estProb.push_back(0.572726279);
//   pyRes->litter.estProb.push_back(0.553298381);
//   pyRes->litter.estProb.push_back(0.543802864);
//   pyRes->litter.estProb.push_back(0.543802864);
//   pyRes->litter.estProb.push_back(0.534476942);
//   pyRes->litter.estProb.push_back(0.534476942);
//   pyRes->litter.estProb.push_back(0.52533783);
//   pyRes->litter.estProb.push_back(0.516402011);
//   pyRes->litter.estProb.push_back(0.516402011);
//
//   pyRes->litter.litterSize.push_back(9);
//   pyRes->litter.litterSize.push_back(9);
//   pyRes->litter.litterSize.push_back(10);
//   pyRes->litter.litterSize.push_back(10);
//   pyRes->litter.litterSize.push_back(11);
//   pyRes->litter.litterSize.push_back(13);
//   pyRes->litter.litterSize.push_back(14);
//   pyRes->litter.litterSize.push_back(14);
//   pyRes->litter.litterSize.push_back(15);
//   pyRes->litter.litterSize.push_back(16);
//   pyRes->litter.litterSize.push_back(9);
//   pyRes->litter.litterSize.push_back(9);
//   pyRes->litter.litterSize.push_back(10);
//   pyRes->litter.litterSize.push_back(10);
//   pyRes->litter.litterSize.push_back(11);
//   pyRes->litter.litterSize.push_back(12);
//   pyRes->litter.litterSize.push_back(13);
//   pyRes->litter.litterSize.push_back(14);
//   pyRes->litter.litterSize.push_back(14);
//   pyRes->litter.litterSize.push_back(14);
//   pyRes->litter.litterSize.push_back(7);
//   pyRes->litter.litterSize.push_back(10);
//   pyRes->litter.litterSize.push_back(10);
//   pyRes->litter.litterSize.push_back(11);
//   pyRes->litter.litterSize.push_back(11);
//   pyRes->litter.litterSize.push_back(11);
//   pyRes->litter.litterSize.push_back(11);
//   pyRes->litter.litterSize.push_back(14);
//   pyRes->litter.litterSize.push_back(14);
//   pyRes->litter.litterSize.push_back(15);
//   pyRes->litter.litterSize.push_back(8);
//   pyRes->litter.litterSize.push_back(10);
//   pyRes->litter.litterSize.push_back(11);
//   pyRes->litter.litterSize.push_back(11);
//   pyRes->litter.litterSize.push_back(12);
//   pyRes->litter.litterSize.push_back(12);
//   pyRes->litter.litterSize.push_back(13);
//   pyRes->litter.litterSize.push_back(14);
//   pyRes->litter.litterSize.push_back(14);   
//
//   pyRes->litter.expected.push_back(1.148272326);
//   pyRes->litter.expected.push_back(1.148272326);
//   pyRes->litter.expected.push_back(1.323471805);
//   pyRes->litter.expected.push_back(1.323471805);
//   pyRes->litter.expected.push_back(1.508194016);
//   pyRes->litter.expected.push_back(1.906206638);
//   pyRes->litter.expected.push_back(2.119497049);
//   pyRes->litter.expected.push_back(2.119497049);
//   pyRes->litter.expected.push_back(2.342310192);
//   pyRes->litter.expected.push_back(2.574646069);
//   pyRes->litter.expected.push_back(2.713743308);
//   pyRes->litter.expected.push_back(2.713743308);
//   pyRes->litter.expected.push_back(2.977817749);
//   pyRes->litter.expected.push_back(2.977817749);
//   pyRes->litter.expected.push_back(3.238103583);
//   pyRes->litter.expected.push_back(3.495536119);
//   pyRes->litter.expected.push_back(3.751018618);
//   pyRes->litter.expected.push_back(4.005420479);
//   pyRes->litter.expected.push_back(4.005420479);
//   pyRes->litter.expected.push_back(4.005420479);
//   pyRes->litter.expected.push_back(3.033741255);
//   pyRes->litter.expected.push_back(4.102322662);
//   pyRes->litter.expected.push_back(4.102322662);
//   pyRes->litter.expected.push_back(4.434656714);
//   pyRes->litter.expected.push_back(4.434656714);
//   pyRes->litter.expected.push_back(4.434656714);
//   pyRes->litter.expected.push_back(4.434656714);
//   pyRes->litter.expected.push_back(5.374621982);
//   pyRes->litter.expected.push_back(5.374621982);
//   pyRes->litter.expected.push_back(5.672427209);
//   pyRes->litter.expected.push_back(4.581810233);
//   pyRes->litter.expected.push_back(5.532983808);
//   pyRes->litter.expected.push_back(5.981831502);
//   pyRes->litter.expected.push_back(5.981831502);
//   pyRes->litter.expected.push_back(6.413723301);
//   pyRes->litter.expected.push_back(6.413723301);
//   pyRes->litter.expected.push_back(6.829391788);
//   pyRes->litter.expected.push_back(7.229628148);
//   pyRes->litter.expected.push_back(7.229628148);
//
//   pyRes->litter.observed.push_back(0);
//   pyRes->litter.observed.push_back(1);
//   pyRes->litter.observed.push_back(1);
//   pyRes->litter.observed.push_back(2);
//   pyRes->litter.observed.push_back(2);
//   pyRes->litter.observed.push_back(3);
//   pyRes->litter.observed.push_back(3);
//   pyRes->litter.observed.push_back(2);
//   pyRes->litter.observed.push_back(2);
//   pyRes->litter.observed.push_back(1);
//   pyRes->litter.observed.push_back(2);
//   pyRes->litter.observed.push_back(5);
//   pyRes->litter.observed.push_back(2);
//   pyRes->litter.observed.push_back(1);
//   pyRes->litter.observed.push_back(4);
//   pyRes->litter.observed.push_back(3);
//   pyRes->litter.observed.push_back(6);
//   pyRes->litter.observed.push_back(6);
//   pyRes->litter.observed.push_back(4);
//   pyRes->litter.observed.push_back(3);
//   pyRes->litter.observed.push_back(2);
//   pyRes->litter.observed.push_back(5);
//   pyRes->litter.observed.push_back(5);
//   pyRes->litter.observed.push_back(4);
//   pyRes->litter.observed.push_back(4);
//   pyRes->litter.observed.push_back(5);
//   pyRes->litter.observed.push_back(4);
//   pyRes->litter.observed.push_back(4);
//   pyRes->litter.observed.push_back(5);
//   pyRes->litter.observed.push_back(6);
//   pyRes->litter.observed.push_back(5);
//   pyRes->litter.observed.push_back(4);
//   pyRes->litter.observed.push_back(6);
//   pyRes->litter.observed.push_back(6);
//   pyRes->litter.observed.push_back(8);
//   pyRes->litter.observed.push_back(8);
//   pyRes->litter.observed.push_back(7);
//   pyRes->litter.observed.push_back(6);
//   pyRes->litter.observed.push_back(6);
//
//   pyRes->litter.SR.push_back(-1.147257987);
//   pyRes->litter.SR.push_back(-0.148141348);
//   pyRes->litter.SR.push_back(-0.301860366);
//   pyRes->litter.SR.push_back(0.631328744);
//   pyRes->litter.SR.push_back(0.431109028);
//   pyRes->litter.SR.push_back(0.857594396);
//   pyRes->litter.SR.push_back(0.65653973);
//   pyRes->litter.SR.push_back(-0.089101984);
//   pyRes->litter.SR.push_back(-0.243481535);
//   pyRes->litter.SR.push_back(-1.071325161);
//   pyRes->litter.SR.push_back(-0.518421334);
//   pyRes->litter.SR.push_back(1.660602955);
//   pyRes->litter.SR.push_back(-0.676196332);
//   pyRes->litter.SR.push_back(-1.367732493);
//   pyRes->litter.SR.push_back(0.504037656);
//   pyRes->litter.SR.push_back(-0.314836859);
//   pyRes->litter.SR.push_back(1.376689417);
//   pyRes->litter.SR.push_back(1.179530167);
//   pyRes->litter.SR.push_back(-0.003205497);
//   pyRes->litter.SR.push_back(-0.594573329);
//   pyRes->litter.SR.push_back(-0.788462565);
//   pyRes->litter.SR.push_back(0.577118305);
//   pyRes->litter.SR.push_back(0.577118305);
//   pyRes->litter.SR.push_back(-0.267167737);
//   pyRes->litter.SR.push_back(-0.267167737);
//   pyRes->litter.SR.push_back(0.347496039);
//   pyRes->litter.SR.push_back(-0.267167737);
//   pyRes->litter.SR.push_back(-0.755412682);
//   pyRes->litter.SR.push_back(-0.205870559);
//   pyRes->litter.SR.push_back(0.174415333);
//   pyRes->litter.SR.push_back(0.298883377);
//   pyRes->litter.SR.push_back(-0.975099884);
//   pyRes->litter.SR.push_back(0.010998302);
//   pyRes->litter.SR.push_back(0.010998302);
//   pyRes->litter.SR.push_back(0.918022312);
//   pyRes->litter.SR.push_back(0.918022312);
//   pyRes->litter.SR.push_back(0.094758157);
//   pyRes->litter.SR.push_back(-0.65761782);
//   pyRes->litter.SR.push_back(-0.65761782);
//
//   pyRes->reduced.dose.push_back(0);
//   pyRes->reduced.dose.push_back(25);
//   pyRes->reduced.dose.push_back(50);
//   pyRes->reduced.dose.push_back(100);
//   pyRes->reduced.propAffect.push_back(0.14049587);
//   pyRes->reduced.propAffect.push_back(0.31034482);
//   pyRes->reduced.propAffect.push_back(0.38596491);
//   pyRes->reduced.propAffect.push_back(0.53333333);
//   pyRes->reduced.lowerConf.push_back(0.10100315);
//   pyRes->reduced.lowerConf.push_back(0.23266263);
//   pyRes->reduced.lowerConf.push_back(0.34156249);
//   pyRes->reduced.lowerConf.push_back(0.46300766);
//   pyRes->reduced.upperConf.push_back(0.19144442);
//   pyRes->reduced.upperConf.push_back(0.39977302);
//   pyRes->reduced.upperConf.push_back(0.43230026);
//   pyRes->reduced.upperConf.push_back(0.60240216);

}


double opt_nlogistic(std::vector<double> &p, const std::vector<double> &Ls, const std::vector<double> &Xi, const std::vector<int> &Xg, const std::vector<double> &Yp, const std::vector<double> &Yn, double smax, double smin, bool isRestricted, int pass){

   std::cout<<"Inside opt_nlogistic"<<std::endl;

   bool fail = true;
   double minf;
   //start with hardcoded parameter limits
   //will move these to header
   double slopeUpperBound = 18.0;
   int nparm = p.size();

   std::vector<double> lb(nparm);
   std::vector<double> ub(nparm);
   //alpha
   lb[0] = 0.0;
   ub[0] = 1.0;
   //beta
   lb[1] = -1.0*DBL_MAX;
   ub[1] = DBL_MAX;
   //theta1
   lb[2] = 0.0;
   if (pass == 1 || pass == 2){
     ub[2] = 0.0;
   } else {
     ub[2] = 1.0;
   }
   //theta2
   if (pass == 1 || pass == 2){
     lb[3] = 0.0;
     ub[3] = 0.0;
   } else {
     lb[3] = -1.0*DBL_MAX;
     ub[3] = DBL_MAX;
   }
   //rho
   if (isRestricted){
     lb[4] = 1.0;
   } else {
     lb[4] = 0.0;
   }
   ub[4] = slopeUpperBound;
   //phi(s)
   for (int i=5; i<nparm; i++){
     lb[i] = 0.0;
     if (pass == 1){
       ub[i] = 0.0;
     } else {
       ub[i] = DBL_MAX;
     }
   }
   
   //Description of optimization problem
   //minimize -log-likelihood
   //inequality constraints
   //	alpha + theta1*rij >=0
   //	alpha + theta1*rij < 1


   nlopt::opt opt(nlopt::LD_SLSQP, nparm);
   //BMD opt
   
   struct nestedObjData objData;
   objData.Ls = Ls;
   objData.Xi = Xi;
   objData.Xg = Xg;
   objData.Yp = Yp;
   objData.Yn = Yn;
   objData.smax = smax;
   objData.smin = smin; 

   opt.set_min_objective(objfunc_nlogistic_ll, &objData);
//   opt.add_inequality_constraint(nlogistic_ineq1, &ineq1, 1e-8);
//   opt.add_inequality_constraint(nlogistic_ineq2, &ineq1, 1e-8);
   opt.set_xtol_rel(1e-8);
   //opt.set_xtol_abs(1e-8);
   opt.set_maxeval(1000000);
   opt.set_lower_bounds(lb);
   opt.set_upper_bounds(ub);

   nlopt::result result = nlopt::FAILURE;

   while (fail){
     try{
       result = opt.optimize(p, minf);
       fail = false;
     } catch (std::exception &e){
       std::cout << "nlogistic opt failed: " << e.what() << std::endl;
     }

   }
  
   std::cout<<"minf = "<<minf<<std::endl; 
   std::cout<<"parm after opt"<<std::endl;
   for(int i=0; i<nparm; i++){
     std::cout<<"i:"<<i<<", p:"<<p[i]<<std::endl;
   } 
   return result;
}


void probability_inrange(double *ex){
  if (*ex < 1.0e-7) *ex = 1.0e-7;
  if (*ex > 0.9999999) *ex = 0.9999999;
}


double objfunc_nlogistic_ll(const std::vector<double> &p, std::vector<double> &grad, void *data){

  nestedObjData *d = reinterpret_cast<nestedObjData*>(data);
  //from data struct
  std::vector<double> Ls = d->Ls;
  std::vector<double> Xi = d->Xi;
  std::vector<int> Xg = d->Xg;
  std::vector<double> Yp = d->Yp;
  std::vector<double> Yn = d->Yn;
  double smax = d->smax;
  double smin = d->smin;

  if (!grad.empty()){
    Nlogist_g(p, Ls, Xi, Xg, Yp, Yn, smax, smin, grad);
  }

  //negative log-likelihood calc
  double loglike = Nlogist_lk(p, Ls, Xi, Xg, Yp, Yn, smax, smin);

  return loglike;

}

// Used to comput the log-likelihood for Nlogistic model
double Nlogist_lk(std::vector<double> p, std::vector<double> Ls, std::vector<double> Xi, std::vector<int> Xg, std::vector<double> Yp, std::vector<double> Yn, double smax, double smin){

   int Nobs = Xi.size();
   int nparms = p.size();
   std::vector<double> probs(Nobs);
   std::vector<std::vector<double>> gradij(Nobs, std::vector<double> (5));  


   bool compgrad = false;
   bool isBMDL = false;
//   std::cout<<"b4 Nlogist_probs"<<std::endl;
   Nlogist_probs(probs, p, compgrad, gradij, isBMDL, smax, smin, Ls, Xi);
//   std::cout<<"after Nlogist_probs"<<std::endl;
//   for(int i=0; i<Nobs; i++){
//     std::cout<<"i:"<<i<<", probs:"<<probs[i]<<std::endl;
//   } 
   
   double tm, tm1, tm2, tm3;
   int plus5, j;
   double xlk = 0.0;   
   for (int i=0; i<Nobs; i++){
     tm1 = 0.0;
     tm2 = 0.0;
     tm3 = 0.0;
     plus5 = 5 + Xg[i];
     j = (int) Yp[i];
     //if (probs[i] == 0.0 && j > 0){ //this should be changed for better floating point comparison
     if (probs[i] < CloseToZero && j > 0){ 
        tm1 -=40.0;
     } else {
       for (int k=1; k<=j; k++){
          tm = probs[i] + (k-1)*p[plus5];
	  tm1 += log(tm);
       }
     }
     j = (int) Yn[i];
     if (probs[i] >=1.0 && j > 0){
       tm2 -=40.0;
     } else {
       for (int k=1; k<=j; k++){
         tm = 1.0 - probs[i] + (k-1)*p[plus5];
	 tm2 += log(tm);
       }
     }
     j = (int) (Yn[i] + Yp[i]);
     for (int k=1; k<=j; k++){
       tm = 1.0 + (k-1)*p[plus5];
       tm3 += log(tm);
     }
     xlk += (tm1 + tm2 - tm3);
   }
   std::cout<<"returning:"<<(-1.0*xlk)<<" from Nlogist_lk"<<std::endl;
  
   //returns negative log likelihood
   return -1.0*xlk;
}


//Used to compute the gradients for Nlogist_model.  Parameters are in "internal" transformed form.
double Nlogist_g(std::vector<double> p, std::vector<double> Ls, std::vector<double> Xi, std::vector<int> Xg, std::vector<double> Yp, std::vector<double> Yn, double smax, double smin, std::vector<double> &g){



  int Nobs = Xi.size();
  int nparm = p.size();
  std::vector<double> probs(Nobs);
  std::vector<std::vector<double>> gradij(Nobs, std::vector<double> (5));
  double ex, tm, tm1, tm2, tm3, tm1a, tm2a, tm3a, tm12;
  int plus5, j;
  std::vector<double> tmp_g(nparm);
  std::vector<double> dd(nparm);


  bool compgrad = true;
  bool isBMDL = false;
  std::cout<<"b4 Nlogist_probs"<<std::endl;
  Nlogist_probs(probs, p, compgrad, gradij, isBMDL, smax, smin, Ls, Xi);


  //initial tmp g[j]'s
  for (int i=0; i<nparm; i++){
    tmp_g[i] = 0.0;
  }
  for (int i=0; i<Nobs; i++){
    ex = probs[i];

    //compute first partial derivatives
    tm1 = tm2 = tm3 = 0.0;
    tm1a = tm2a = tm3a = 0.0;
    for (int j=5; j<nparm; j++){
      dd[j] = 0.0;
    }    
    plus5 = 5 + Xg[i];
    j = (int) Yp[i];
    if (ex > 0.0){
      for (int k=1; k<=j; k++){
        tm = ex + (k-1)*p[plus5];
	tm1 += 1.0/tm;
	tm1a += (1.0/tm)*(k-1);
      }
    }
    j = (int) Yn[i];
    if (ex < 1.0){
      for (int k=1; k<=j; k++){
        tm = 1.0-ex+(k-1)*p[plus5];
	tm2 += 1.0/tm;
	tm2a += (1.0/tm)*(k-1);
      }
    }
    j = (int) (Yn[i] + Yp[i]);
    for (int k=1; k<=j; k++){
      tm = 1.0+(k-1)*p[plus5];
      if (tm == 0.0) tm = 0.000001;
      tm3 += 1.0/tm;
      tm3a += (1.0/tm)*(k-1);
    }
    tm12 = (tm1-tm2);
    for (int k=0; k<5; k++){
      dd[k] = gradij[i][k]*tm12;
    }
    for (int k=0; k<nparm; k++){
      tmp_g[k] -= dd[k];
    }
  }
  //end of 1st partial derivative
  g = tmp_g;
 
  return 0;
}

void Nlogist_probs(std::vector<double> &probs, const std::vector<double> &p, bool compgrad, std::vector<std::vector<double>> &gradij, bool isBMDL, double smax, double smin, const std::vector<double> &Ls, const std::vector<double> &Xi){

  double spij, smij, ex, ex1, ex2, ex3, dd2;

  double sdif = smax - smin;
  int nparms = p.size();
  int Nobs = Xi.size();

  std::vector<double> pint (nparms);
  for (int i=0; i<nparms; i++){
    pint[i] = p[i];
  }


  if (isBMDL){
	  ////FILL THIS IN
  }

  for (int i=0; i<Nobs; i++){
    spij = (smax - Ls[i])/sdif;
    smij = (Ls[i] - smin)/sdif;

    //enforce alpha+Theta1*Sij >= 0
    //We currently don't allow the user to specify Theta parameters
    //will need to possibly address this for python code
    //if (Spec[0] == Spec[2]){
    //   ex = spij * pint[0] + smij * pint[2];
    //} else {
       ex = pint[0] + Ls[i] * pint[2];
    //}
    ex2 = 1.0 - ex;
    ex1 = 0.0;
    if (Xi[i] > 0){
       ex1 = exp(-(pint[1] + pint[3] * Ls[i] + pint[4] * log(Xi[i])));
       ex = ex + ex2/(1+ex1);
    }
    probability_inrange(&ex);
    probs[i] = ex;    

    if (compgrad){
       ex3 = ex1 / ((1.0+ex1)*(1.0+ex1));
       dd2 = ex2*ex3;

       if (isBMDL){
	 ////FILL THIS IN
       } else {
         //enforce alpha+Theta1*Sij >= 0
         //We currently don't allow specifying Theta parameters
         //will need to possibly address this for python code
	 //if (Spec[0] == Spec[2]){
         //   if (Xi[i] > 0){
	 //     gradij]i][0] = spij * (1.0 - 1.0/(1.0+ex1));
	 //     gradij[i][2] = smij * (1.0 - 1.0/(1.0+ex1));
	 //   } else {
	 //     //Case where Xi[i] == 0
         //     gradij[i][0] = spij;
	 //     gradij[i][2] = smij;
	 //   }
	 //} else {
            if (Xi[i] > 0){
	      gradij[i][0] = 1.0 - 1.0/(1.0+ex1);
	      gradij[i][2] = Ls[i] * gradij[i][0];
	    } else {
	      //Case where Xi[i] == 0
	      gradij[i][0] = 1.0;
	      gradij[i][2] = Ls[i];
	    }
	 //} //end of derivates that depend on the transformation of alpha and theta1
	 if (Xi[i] > 0){
           gradij[i][1] = dd2;
	   gradij[i][3] = dd2 * Ls[i];
	   gradij[i][4] = dd2 * log(Xi[i]);
	 } else {
	   gradij[i][1] = 0.0;
	   gradij[i][3] = 0.0;
	   gradij[i][4] = 0.0;
	 }
       }
    }

  }

};
