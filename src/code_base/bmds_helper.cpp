#include <stdio.h>
#include <math.h>
//#include <cmath>
#include "bmds_helper.h"
#include "analysis_of_deviance.h"


int checkForBoundedParms(int nparms, double *parms, double *lowerBound, double *upperBound, struct BMDS_results *BMDSres ){
   // First find number of bounded parms
   int bounded = 0;
//   std::cout << "checkForBoundedParms with tol:"<<BMDS_EPS<<"\n";
   for (int i=0; i<nparms; i++){
      //5*i+4 is location of min in prior array
      //5*i+5 is location of max in prior array 
//      std::cout<<"i:"<<i<<", parms[i]: "<<parms[i]<<"\n";
//      std::cout<<"lowerBound: " << lowerBound[i] << std::endl;
//      std::cout<<"upperBound: " << upperBound[i] << std::endl;
//      std::cout<<"lower bound check: "<<fabs(parms[i]-lowerBound[i])<<"\n";
//      std::cout<<"upper bound check: "<<fabs(parms[i]-upperBound[i])<<"\n";
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

}


void collect_dichoMA_bmd_values(struct dichotomousMA_analysis *anal, struct dichotomousMA_result *res, struct BMDSMA_results *BMDSres){

  int distSize = res->dist_numE*2;

  double* quant = (double*)malloc(distSize / 2 * sizeof(double));
  double* val = (double*)malloc(distSize / 2 * sizeof(double));

  for (int i = 0; i < distSize/2; i++){
    val[i] = res->bmd_dist[i];
  }
  for (int i = distSize/2; i < distSize; i++){
    quant[i-distSize/2] = res->bmd_dist[i];
  }

//  calculate MA quantiles
  BMDSres->BMD_MA = findQuantileVals(quant, val, distSize/2, 0.50);
  BMDSres->BMDL_MA = findQuantileVals(quant, val, distSize/2, 0.05);
  BMDSres->BMDU_MA = findQuantileVals(quant, val, distSize/2, 0.95);

// calculate individual model quantiles
  for (int j=0; j<anal->nmodels; j++){
      for (int i = 0; i < distSize/2; i++){
        val[i] = res->models[j]->bmd_dist[i];
      }
      for (int i = distSize/2; i < distSize; i++){
        quant[i-distSize/2] = res->models[j]->bmd_dist[i];
      }
      BMDSres->BMD[j] = findQuantileVals(quant, val, distSize/2, 0.50);
      BMDSres->BMDL[j] = findQuantileVals(quant, val, distSize/2, 0.05);
      BMDSres->BMDU[j] = findQuantileVals(quant, val, distSize/2, 0.95);
  }
}

void collect_cont_bmd_values(struct continuous_analysis *anal, struct continuous_model_result *res, struct BMDS_results *BMDSres){

  int distSize = res->dist_numE*2;

  double* quant = (double*)malloc(distSize / 2 * sizeof(double));
  double* val = (double*)malloc(distSize / 2 * sizeof(double));

  for (int i = 0; i < distSize/2; i++){
    val[i] = res->bmd_dist[i];
  }
  for (int i = distSize/2; i < distSize; i++){
    quant[i-distSize/2] = res->bmd_dist[i];
  }

  calcContAIC(anal, res, BMDSres);
  BMDSres->BMD = findQuantileVals(quant, val, distSize/2, 0.50);
  BMDSres->BMDL = findQuantileVals(quant, val, distSize/2, 0.05);
  BMDSres->BMDU = findQuantileVals(quant, val, distSize/2, 0.95);

  
}


void BMDS_ENTRY_API __stdcall runBMDSDichoAnalysis(struct dichotomous_analysis *anal, struct dichotomous_model_result *res, struct dichotomous_GOF *gof, struct BMDS_results *bmdsRes, struct dicho_AOD *bmdsAOD){

  bmdsRes->validResult = false;

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
    gof->expected[i] = gofRes.expected[i];
	gof->residual[i] = gofRes.residual[i];
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
    gof->ebLower[i] = (eb1 - 1 - eb2)/ ebDenom; 
    gof->ebUpper[i] = (eb1 + 1 + eb2)/ ebDenom;
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
    bmdsRes->stdErr[i] = BMDS_MISSING;
    bmdsRes->lowerConf[i] = BMDS_MISSING;
    bmdsRes->upperConf[i] = BMDS_MISSING;
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
 
}


void BMDS_ENTRY_API __stdcall runBMDSContAnalysis(struct continuous_analysis *anal, struct continuous_model_result *res, struct BMDS_results *bmdsRes, struct continuous_AOD *aod, struct continuous_GOF *gof, bool *detectAdvDir, bool *restricted){

  bmdsRes->validResult = false;
  anal->transform_dose = false;

  //debug output
  std::cout.precision(17);
  std::cout<<"detectAdvDir:" << (*detectAdvDir ? "true":"false")<<std::endl;
  std::cout<<"restricted:" << (*restricted ? "true":"false")<<std::endl;
  std::cout<<"continuous_analysis struct:" << std::endl;
  std::cout<<"---------------------------"<<std::endl;
  std::cout<<"cont_model:"<<anal->model<<std::endl;
  std::cout<<"n:"<<anal->n<<std::endl;
  std::cout<<"suff_stat:" << (anal->suff_stat ? "true":"false")<<std::endl;
  std::cout<<"BMD_type:" << anal->BMD_type<<std::endl;
  std::cout<<"isIncreasing:"<<(anal->isIncreasing ? "true":"false")<<std::endl;
  std::cout<<"BMR:"<<anal->BMR<<std::endl;
  std::cout<<"disttype:"<<anal->disttype<<std::endl;
  std::cout<<"alpha:"<<anal->alpha<<std::endl;
  std::cout<<"degree:"<<anal->degree<<std::endl;
  std::cout<<"parms:"<<anal->parms<<std::endl;
  std::cout<<"prior_cols:"<<anal->prior_cols<<std::endl;
  std::cout<<"transform_dose:"<<anal->transform_dose<<std::endl;
  std::cout<<"Data:"<<std::endl;
  std::cout<<"Dose, N, Mean, Std. Dev"<<std::endl;
  for (int i=0; i<anal->n; i++){
    std::cout<< anal->doses[i]<<", "<<anal->n_group[i]<<", "<<anal->Y[i]<<",  "<<anal->sd[i]<<std::endl;
  }
  std::cout<<"Parameter settings/Priors"<<std::endl;
  for (int i=0; i<anal->prior_cols*anal->parms; i++){
    std::cout<<anal->prior[i]<<std::endl;
  }


  //if (anal->model == cont_model::polynomial && anal->disttype == distribution::log_normal){
  if (anal->model != cont_model::exp_3 && anal->model != cont_model::exp_5){
    if(anal->disttype == distribution::log_normal){
      std::cout << "lognormal distribution is only compatible with exponential models\n";
      return; 
    }
  }

  if (*detectAdvDir){
    std::cout << "IN detectAdvDir\n";
    determineAdvDir(anal);

    int ind;
    if (*restricted) {
      std::cout << "IN restricted\n";
      switch(anal->model){
//        case cont_model::exp_3:
//        case cont_model::exp_5:
//          if (anal->prior[0] == 0){   //checks if frequentist model
//            if(anal->isIncreasing){
//              if (anal->disttype == distribution::normal_ncv){
//                anal->prior[20] = 0.0;  //c min
//              } else {
//                anal->prior[17] = 0.0;  //c min
//              }
//            } else {
//              if (anal->disttype == distribution::normal_ncv){
//                anal->prior[26] = 0.0;  //c max
//              } else {
//                anal->prior[22] = 0.0;  //c max
//              }
//            }
//          }
//          break;
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
//            //set ind for alpha parm
//            if(anal->isIncreasing){
//              ind = 4*numRows-1;  
//            } else {
//              ind = 5*numRows-1;
//            }
//            anal->prior[ind] = 0.0;  //alpha min or max
          }
          break;
//        case cont_model::hill:
//          if (anal->isIncreasing){
//             if (anal->disttype == distribution::normal_ncv){
//                anal->prior[20] = 0.0;
//             } else {
//                anal->prior[17] = 0.0;
//             }
//          } else {
//             if (anal->disttype == distribution::normal_ncv){
//                anal->prior[26] = 0.0;
//             } else {
//                anal->prior[22] = 0.0;
//             }
//          }
//          break; 
      }  //end switch
    } //end if restricted
  } //end if detectAdvDir


  std::cout<<"Parameter settings/Priors POST changes"<<std::endl;
  for (int i=0; i<anal->prior_cols*anal->parms; i++){
    std::cout<<anal->prior[i]<<std::endl;
  }

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

  //tmp fix for exp_3 parameter shift
//  if (anal->model == cont_model::exp_3){
//    struct continuous_model_result tmpRes;
//    double* tmpParms = (double*)malloc(res->nparms * sizeof(double));
//    double* tmpCov = (double*)malloc(res->nparms*res->nparms * sizeof(double));
//    double* tmpBmd_dist = (double*)malloc(res->dist_numE*2 * sizeof(double));
//
//    tmpRes.model = res->model;
//    tmpRes.dist = res->dist;
//    tmpRes.nparms = res->nparms;
//    tmpRes.max = res->max;
//    tmpRes.dist_numE = res->dist_numE;
//    tmpRes.model_df = res->model_df;
//    tmpRes.total_df = res->total_df;
//    tmpRes.bmd = res->bmd;
//
//    tmpRes.parms = tmpParms;
//    for(int i=0; i<res->nparms; i++){
//      tmpRes.parms[i] = res->parms[i]; 
//    }
//    //shift all parms after b to higher pos
//    for(int i=res->nparms-1; i>1; i--){
//      tmpRes.parms[i] = tmpRes.parms[i-1];
//    }
//
//    tmpRes.cov = tmpCov;
//    for(int i=0; i<res->nparms*res->nparms; i++){
//       tmpRes.cov[i] = res->cov[i];
//    }
//    tmpRes.bmd_dist = tmpBmd_dist;
//    for(int i=0; i<res->dist_numE*2; i++){
//      tmpRes.bmd_dist[i] = res->bmd_dist[i];
//    }
//
//
//    continuous_expectation(&GOFanal, &tmpRes, &GOFres);
//  } else {
    continuous_expectation(&GOFanal, res, &GOFres);
//  }
 

  for (int i=0; i<GOFanal.n; i++){
	gof->dose[i] = GOFanal.doses[i];
    gof->size[i] = GOFanal.n_group[i];
    gof->estMean[i] = GOFres.expected[i];
    gof->obsMean[i] = GOFanal.Y[i];
    gof->estSD[i] = GOFres.sd[i];
    gof->obsSD[i] = GOFanal.sd[i];
    gof->res[i] = sqrt(gof->size[i])*(gof->obsMean[i] - gof->estMean[i]) / gof->estSD[i];
  }
  if (anal->disttype == distribution::log_normal){
    for (int i=0; i<GOFanal.n; i++){
      gof->calcMean[i] = exp(log(GOFanal.Y[i]) - log(1 + pow(GOFanal.sd[i] / GOFanal.Y[i], 2.0)) / 2);
      gof->calcSD[i] = exp(sqrt(log(1.0 + pow(GOFanal.sd[i]/GOFanal.Y[i], 2.0))));
    }
    for (int i=0; i<GOFanal.n; i++){
      gof->estMean[i] = exp(gof->estMean[i]+pow(exp(res->parms[res->nparms-1]),2)/2);
      gof->res[i] = sqrt(gof->size[i])*(gof->obsMean[i] - gof->estMean[i]) / gof->estSD[i];
    }
  } else {
	for (int i=0; i<GOFanal.n; i++){
	  gof->calcMean[i] = GOFanal.Y[i];
          gof->calcSD[i] = GOFanal.sd[i];
	}
  }

  for (int i=0; i<GOFanal.n; i++){
    gof->ebLower[i] = gof->calcMean[i] + gsl_cdf_tdist_Pinv(0.025, gof->n - 1) * (gof->obsSD[i]/sqrt(gof->size[i]));
    gof->ebUpper[i] = gof->calcMean[i] + gsl_cdf_tdist_Pinv(0.975, gof->n - 1) * (gof->obsSD[i]/sqrt(gof->size[i]));
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
  
  calc_contAOD(anal, &GOFanal, res, bmdsRes, aod);

  rescale_contParms(anal, res->parms); 

  for (int i=0; i< anal->parms; i++){
    bmdsRes->stdErr[i] = BMDS_MISSING;
    bmdsRes->lowerConf[i] = BMDS_MISSING;
    bmdsRes->upperConf[i] = BMDS_MISSING;
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
  
//  struct dichotomous_aod aod;
//  double A1 = BMDS_MISSING;
//  int N1 = BMDS_MISSING;
//  double A2 = BMDS_MISSING;
//  int N2 = BMDS_MISSING;
//  aod.A1 = A1;
//  aod.N1 = N1;
//  aod.A2 = A2;
//  aod.N2 = N2;

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

  bmdsAOD->devRed = dev = 2* (bmdsAOD->fittedLL - bmdsAOD->redLL);
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
      std::cout<<"lognormal addConst calcs"<<std::endl;
      for (int i=0; i<CA->n;i++){
         std::cout<<"i:"<<i<<", Y:"<<CA->Y[i]<<std::endl;
      }
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

  struct testsOfInterest *TOI = aod->TOI;
  double dev;
  int df;

  //Test #1 - A2 vs Reduced - does mean and/or variance differ across dose groups
  //TOI->llRatio[0] = dev = 2 * (aod->LL[1] - aod->LL[4]);
  dev = (aod->LL[1] - aod->LL[4]);
  //handle underflow/negative zero
  if (dev < BMDS_EPS){
     dev = 0.0;
  }
  TOI->llRatio[0] = dev = 2*dev;
  TOI->DF[0] = df = aod->nParms[1] - aod->nParms[4];
  TOI->pVal[0] = (isnan(dev) || dev < 0.0 || df < 0.0) ? -1 : 1.0 - gsl_cdf_chisq_P(dev, df);

  //Test #2 - A1 vs A2 - homogeneity of variance across dose groups
  dev = (aod->LL[1] - aod->LL[0]);
  //handle underflow/negative zero
  if (dev < BMDS_EPS){
     dev = 0.0;
  }
  TOI->llRatio[1] = dev = 2*dev;
  //TOI->llRatio[1] = dev = 2 * (aod->LL[1] - aod->LL[0]);
  TOI->DF[1] = df = aod->nParms[1] - aod->nParms[0];
  TOI->pVal[1] = (isnan(dev) || dev < 0.0 || df < 0.0) ? -1 : 1.0 - gsl_cdf_chisq_P(dev, df);

  //Test #3 - A2 vs A3 - Does the model describe variances adequately
  dev = (aod->LL[1] - aod->LL[2]);
  //handle underflow/negative zero
  if (dev < BMDS_EPS){
     dev = 0.0;
  }
  //TOI->llRatio[2] = dev = 2 * (aod->LL[1] - aod->LL[2]);
  TOI->llRatio[2] = dev = 2*dev;
  TOI->DF[2] = df = aod->nParms[1] - aod->nParms[2];
  TOI->pVal[2] = (isnan(dev) || dev < 0.0 || df < 0.0) ? -1 : 1.0 - gsl_cdf_chisq_P(dev, df);

  //Test #4 - A3 vs Fitted - does the fitted model describe the obs data adequately 
  dev = (aod->LL[2] - aod->LL[3]);
  //handle underflow/negative zero
  if (dev < BMDS_EPS){
     dev = 0.0;
  }
  //TOI->llRatio[3] = dev = 2 * (aod->LL[2] - aod->LL[3]);
  TOI->llRatio[3] = dev = 2*dev;
  TOI->DF[3] = df = aod->nParms[2] - aod->nParms[3];
  TOI->pVal[3] = (isnan(dev) || dev < 0.0 || df < 0.0) ? -1 : 1.0 - gsl_cdf_chisq_P(dev, df);

  //DOF check for test 4
  if (TOI->DF[3] <= 0) {
     TOI->pVal[3] = BMDS_MISSING;
  }
}

void determineAdvDir(struct continuous_analysis *CA){

  int n_rows;

  struct continuous_analysis CAnew;
  //allocate memory for individual size.  Will not use entire array for summary data.
  double* doses = (double*)malloc(CA->n * sizeof(double));
  double* means = (double*)malloc(CA->n * sizeof(double));
  double* n_group = (double*)malloc(CA->n * sizeof(double));
  double* sd = (double*)malloc(CA->n * sizeof(double));
 
  CAnew.doses = doses;
  CAnew.Y = means;
  CAnew.n_group = n_group;
  CAnew.sd = sd;

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
      cleanDouble(&aod->TOI->llRatio[j]);
      cleanDouble(&aod->TOI->DF[j]);
      cleanDouble(&aod->TOI->pVal[j]);
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


void BMDS_ENTRY_API __stdcall version(char * versionStr){
  strcpy(versionStr, BMDS_VERSION);   
}

