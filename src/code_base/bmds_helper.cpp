#ifdef WIN32
	#include "pch.h"
#else
	#include "stdafx.h"
#endif
#include <stdio.h>
#include <iomanip>
#include <math.h>
//#include <cmath>
#include <nlopt.hpp>
#include "bmds_helper.h"
#include "analysis_of_deviance.h"

// calendar versioning; see https://peps.python.org/pep-0440/#pre-releases
std::string BMDS_VERSION = "2023.10a1";


double python_dichotomous_model_result::getSRAtDose(double targetDose, std::vector<double> doses){
   std::vector<double> diff;
   double absDiff = DBL_MAX;
   double srVal = BMDS_MISSING;
   if(!bmdsRes.validResult ||  doses.size() != gof.residual.size()){
      return BMDS_MISSING;
   }
   for (int i=0; i<doses.size(); i++){
      diff.push_back(abs(targetDose - doses[i]));
   }
   int minIndex = std::distance(std::begin(diff), std::min_element(std::begin(diff), std::end(diff)));

   return gof.residual[minIndex];
}


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
  BMDSres->BMDL = findQuantileVals(quant, val, distSize/2, anal->alpha);
  BMDSres->BMDU = findQuantileVals(quant, val, distSize/2, 1.0-anal->alpha);

  free(quant);
  free(val);

}


void collect_dichoMA_bmd_values(struct dichotomousMA_analysis *anal, struct dichotomousMA_result *res, struct BMDSMA_results *BMDSres, double alpha){

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
  BMDSres->BMDL_MA = findQuantileVals(quantMA, valMA, distSize/2, alpha);
  BMDSres->BMDU_MA = findQuantileVals(quantMA, valMA, distSize/2, 1.0-alpha);

// calculate individual model quantiles
  for (int j=0; j<anal->nmodels; j++){
      for (int i = 0; i < distSize/2; i++){
        valMA[i] = res->models[j]->bmd_dist[i];
      }
      for (int i = distSize/2; i < distSize; i++){
        quantMA[i-distSize/2] = res->models[j]->bmd_dist[i];
      }
      BMDSres->BMD[j] = findQuantileVals(quantMA, valMA, distSize/2, 0.50);
      BMDSres->BMDL[j] = findQuantileVals(quantMA, valMA, distSize/2, alpha);
      BMDSres->BMDU[j] = findQuantileVals(quantMA, valMA, distSize/2, 1.0-alpha);
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
  BMDSres->BMDL = findQuantileVals(contQuant, contVal, distSize/2, anal->alpha);
  BMDSres->BMDU = findQuantileVals(contQuant, contVal, distSize/2, 1.0-anal->alpha);

  free(contVal);
  free(contQuant);
}

/*****************************************************************
 * DLgama - double log gamma
 *  input:
 *   X - double value
 *  Returns log(z-1)!
 *************************/
double DLgamma(double x){

   double r, az, z, y, dl;
   if (x==0) return 1;
   z = x;
   r = 12.0;
   az = 1.0;
   for (;z-r < 0; z++) az = az * z;
  y = z * z;
  dl=0.91893853320467274-z+(z-0.5)*log(z)+z*
     (0.08333333333333333*y+0.0648322851140734)/
     (y*(y+0.811320754825416)+0.01752021563249146);
  dl = dl - log(az);
  return dl ;  

}

//calculate dichotomous log likelihood constant
double LogLik_Constant(std::vector<double> Y, std::vector<double> n_group){
   
   double X, N, NX, con, xf, nf, nxf, val;
   
   con = 0;
   for (int i=0; i<Y.size(); i++){
      X = Y[i]+1;
      N = n_group[i] + 1;
      NX = n_group[i] - Y[i] + 1;
      xf = exp(DLgamma(X));
      nf = exp(DLgamma(N));
      nxf = exp(DLgamma(NX));
      val = log(nf/(xf*nxf));
      con += val;
   }
   
   return con;
}




/*
 *  ************************************************************************
 *	    		    C math library
 * function ZEROIN - obtain a function zero within the given range
 *
 * Input
 *	double zeroin(ax,bx,tol, f,nparm, parm, gtol)
 *	double ax; 			Root will be sought for within
 *	double bx;  			a range [ax,bx]
 *	double (*f)(nparm, parm, double x, double gtol); Name of the function whose zero
 *					will be sought
 *	double tol;			Acceptable tolerance for the root value.
 *					May be specified as 0.0 to cause
 *					the program to find the root as
 *					accurate as possible
 *     int nparm                        length of parameter vector to f
 *     double parm[]                    vector of parameters to f
 *     double gtol                      additional scaler parameter to f
 *
 * Output
 *	Zeroin returns an estimate for the root with accuracy
 *	4*EPSILON*abs(x) + tol
 *
 * Algorithm
 *	G.Forsythe, M.Malcolm, C.Moler, Computer methods for mathematical
 *	computations. M., Mir, 1980, p.180 of the Russian edition
 *
 *	The function makes use of the bissection procedure combined with
 *	the linear or quadric inverse interpolation.
 *	At every step program operates on three abscissae - a, b, and c.
 *	b - the last and the best approximation to the root
 *	a - the last but one approximation
 *	c - the last but one or even earlier approximation than a that
 *		1) |f(b)| <= |f(c)|
 *		2) f(b) and f(c) have opposite signs, i.e. b and c confine the root
 *	At every step Zeroin selects one of the two new approximations, the
 *	former being obtained by the bissection procedure and the latter
 *	resulting in the interpolation (if a,b, and c are all different
 *	the quadric interpolation is utilized, otherwise the linear one).
 *	If the latter (i.e. obtained by the interpolation) point is 
 *	reasonable (i.e. lies within the current interval [b,c] not being
 *	too close to the boundaries) it is accepted. The bissection result
 *	is used in the other case. Therefore, the range of uncertainty is
 *	ensured to be reduced at least by the factor 1.6
 *
 ************************************************************************
 */

double zeroin(double ax,double bx, double tol,
	      double (*f)(int, double [], double, double), int nparm,
	      double Parms[], double ck)		
     /* ax        Left border | of the range */
     /* bx        Right border | the root is sought*/
     /* f	  Function under investigation */
     /* nparm     number of parameters in Parms */
     /* Parms     vector of parameters to pass to f */
     /* tol       Acceptable tolerance for the root */
     /* gtol      tolerance to pass to f */
//TODO: clean up comments after debugging is complete
{
  double a,b,c;				/* Abscissae, descr. see above	*/
  double fa;				/* f(a)				*/
  double fb;				/* f(b)				*/
  double fc;				/* f(c)				*/

  std::cout<< std::fixed << std::showpoint;
  std::cout << std::setprecision(15); 
//  std::cout<<"inside zeroin"<<std::endl;
//  std::cout<<"ax="<<ax<<", bx="<<bx<<", tol="<<tol<<std::endl; 
//  std::cout<<"nparm="<<nparm<<std::endl;
//  int i;
//  for (i=0;i<nparm;i++){
//	  std::cout<<"i="<<i<<", Parms[i]="<<Parms[i]<<std::endl;
//  }

  a = ax;  b = bx;
//  printf("fa calc\n");
  fa = (*f)(nparm-1, Parms, a, ck);
//  printf("fb calc\n");
  fb = (*f)(nparm-1, Parms, b, ck);
//  printf("fa=%g\n", fa);
//  printf("fb=%g\n", fb);
  c = a;   fc = fa;

  for(;;)		/* Main iteration loop	*/
  {
    double prev_step = b-a;		/* Distance from the last but one*/
					/* to the last approximation	*/
    double tol_act;			/* Actual tolerance		*/
    double p;      			/* Interpolation step is calcu- */
    double q;      			/* lated in the form p/q; divi- */
  					/* sion operations is delayed   */
 					/* until the last moment	*/
    double new_step;      		/* Step at this iteration       */
//    std::cout<<"start of loop"<<std::endl;
//    std::cout<<"a="<<a<<", b="<<b<<", c="<<c<<std::endl;
//    std::cout<<"fa="<<fa<<", fb="<<fb<<std::endl;
    if( fabs(fc) < fabs(fb) )
    {                         		/* Swap data for b to be the 	*/
	a = b;  b = c;  c = a;          /* best approximation		*/
	fa=fb;  fb=fc;  fc=fa;
    }
    tol_act = 2*DBL_EPSILON*fabs(b) + tol/2;
    new_step = (c-b)/2;

//    std::cout<<"tol_act="<<tol_act<<std::endl;
    if( fabs(new_step) <= tol_act || fb == (double)0 ){
//      std::cout<<"returning b:"<<b<<std::endl;
      return b;				/* Acceptable approx. is found	*/
    }

//    std::cout<<"continuing for another loop"<<std::endl;
    			/* Decide if the interpolation can be tried	*/
    if( fabs(prev_step) >= tol_act	/* If prev_step was large enough*/
	&& fabs(fa) > fabs(fb) )	/* and was in true direction,	*/
    {					/* Interpolatiom may be tried	*/
//	std::cout<<"trying interpolation"<<std::endl;
	double t1,cb,t2;
	cb = c-b;
	if( a==c )			/* If we have only two distinct	*/
	{				/* points linear interpolation 	*/
	  t1 = fb/fa;			/* can only be applied		*/
	  p = cb*t1;
	  q = 1.0 - t1;
 	}
	else				/* Quadric inverse interpolation*/
	{
	  q = fa/fc;  t1 = fb/fc;  t2 = fb/fa;
	  p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
	  q = (q-1.0) * (t1-1.0) * (t2-1.0);
	}
	if( p>(double)0 )		/* p was calculated with the op-*/
	  q = -q;			/* posite sign; make p positive	*/
	else				/* and assign possible minus to	*/
	  p = -p;			/* q				*/

	if( p < (0.75*cb*q-fabs(tol_act*q)/2)	/* If b+p/q falls in [b,c]*/
	    && p < fabs(prev_step*q/2) )	/* and isn't too large	*/
	  new_step = p/q;			/* it is accepted	*/
					/* If p/q is too large then the	*/
					/* bissection procedure can 	*/
					/* reduce [b,c] range to more	*/
					/* extent			*/
    }

    if( fabs(new_step) < tol_act )	/* Adjust the step to be not less*/
      {
//	std::cout<<"adjusting step"<<std::endl;
	if( new_step > (double)0 )	/* than tolerance		*/
	  new_step = tol_act;
	else
	  new_step = -tol_act;
      }

    a = b;  fa = fb;			/* Save the previous approx.	*/
    b += new_step;
//	printf("2nd fb calc\n");
    fb = (*f)(nparm-1, Parms, b, ck);	/* Do step to a new approxim.	*/
    if( (fb > 0 && fc > 0) || (fb < 0 && fc < 0) )
    {                 			/* Adjust c for it to have a sign*/
      c = a;  fc = fa;                  /* opposite to that of b	*/
    }
//    std::cout<<"a="<<a<<", b="<<b<<", c="<<c<<std::endl;
//    std::cout<<"fa="<<fa<<", fb="<<fb<<std::endl;
  }

}


/*****************************************************************
 * BMD_func -- used to compute the values of functions BMD_f at
 *            the point x, given the parm p[] and number of parm.
 *            (ck is another parameter).
 *            This routine is called by Binary_root().
 *  input: n is the number of parameters
 *         p[] is a vector of parameters
 *         x is the natural log of a  dose level
 *         ck
 *  output: value of the function
 *****************************************************************/
double BMD_func(int n, double p[], double x, double ck)
{
  double  poly, fx, D;
  int j;
//  printf("\ninside BMD_func\n");
  D = exp(x);
//  printf("D=%g\n",D);
  poly=p[n];
//  printf("poly=%g\n",poly);
  for (j=n-1; j>0; j--) {
//     printf("j=%d, poly b4=%g\n", j, poly);
     poly=poly*D+p[j];
//     printf("poly after=%g\n", poly);
  }
  poly = poly*D;
//  printf("final poly=%g\n",poly);
//  printf("ck = %g\n", ck);
  fx = poly - ck; /* ck = log(1-A) */
  return fx;
}


double getclmt(python_multitumor_analysis *pyAnal, python_multitumor_result *pyRes, double Dose, double target, double maxDose, std::vector<double> xParms, bool isBMDL){
 
   int nT = pyRes->selectedModelIndex.size(); 
   double bmr = pyAnal->BMR;
   std::vector<int> degree;
   for (int j=0; j<nT; j++){
     int modDeg = pyAnal->models[j][pyRes->selectedModelIndex[j]].degree;
     degree.push_back(modDeg);
   } 
   int N = xParms.size() + 1;
   double bmd = Dose;
   std::vector<double> x(N);
   int iOffset = 0;
   x[0] = log(bmd);

   std::vector<std::vector<double>> tmp2;
   //restructure to group like terms together
   int count = 0;
   for (int i=0; i<nT; i++){
      std::vector<double> theParms;
      for (int j=0; j<=degree[i]; j++){
	 theParms.push_back(xParms[count]);
	 count++;
      }
      tmp2.push_back(theParms);
   }

   int maxDegree = *max_element(degree.begin(), degree.end());
   count = 1;
   for (int j=0; j<=maxDegree; j++){
      for (int i=0; i<nT; i++){
	 if (j<tmp2[i].size()){
            x[count] = tmp2[i][j];
	    count++;
	 }
      }
   } 

   //need to round to roughly single-precision for convergence?????
   for (int j=1; j<x.size(); j++){
      x[j] = round_to(x[j], 0.00001);
   }

   std::vector<double> lb(x.size());
   std::vector<double> ub(x.size());

   //calculate the values that will correspond to '0' and 'infinity' in BMD calcs
   double lminbmd = log(DBL_MIN) - log(maxDose);
   double lmaxbmd = log(DBL_MAX) - log(maxDose);

   lb[0] = lminbmd;  //BMD lower limit
   for (int i=1; i<x.size(); i++){
     lb[i] = 0.0; //beta min value
   }

   ub[0] = log(maxDose);
   for (int i=1; i<x.size(); i++){
     ub[i] = 1e4; //beta max value
   }

   
//   //constraint data
   struct msComboEq eq1;
   eq1.bmr = bmr;
   eq1.nT = nT;
   eq1.degree = degree;

   struct msComboInEq ineq1;
   ineq1.nT = nT;
   ineq1.target = target;
   //TODO need to add handling of failed datasets

   for(int i=0; i<pyRes->selectedModelIndex.size(); i++){
     int selIndex = pyRes->selectedModelIndex[i];
     std::vector<double> scaledDose = pyAnal->models[i][selIndex].doses;
     for (int j=0; j<scaledDose.size(); j++){
       scaledDose[j] /= maxDose;
     }
     ineq1.doses.push_back(scaledDose);
     ineq1.Y.push_back(pyAnal->models[i][selIndex].Y);
     ineq1.n_group.push_back(pyAnal->models[i][selIndex].n_group);
   }
   ineq1.nObs = pyAnal->n;
   ineq1.degree = degree;


   nlopt::opt opt(nlopt::LD_SLSQP, x.size());
   if (isBMDL){
      opt.set_min_objective(objfunc_bmdl, NULL);
   } else {
      opt.set_min_objective(objfunc_bmdu, NULL);
   }
   opt.add_equality_constraint(myEqualityConstraint, &eq1, 1e-8);
   opt.add_inequality_constraint(myInequalityConstraint1, &ineq1, 1e-8);
   opt.set_xtol_rel(1e-8);
   opt.set_maxeval(1000000);
   opt.set_lower_bounds(lb);
   opt.set_upper_bounds(ub);

   double minf, val;
   nlopt::result result = nlopt::FAILURE;
   try{
     result = opt.optimize(x, minf);
//     std::cout << "found minimum at f(" << x[0] << ") = " << std::setprecision(10) << minf << std::endl;
   } catch (std::exception &e){
     std::cout << "nlopt failed: " << e.what() << std::endl;
   }

   val = x[0];


   return val;

}


/*****************************************************************
 *  Added by CVL 7/2007 - start block
 * BMDL_combofunc -- returns the lower confidence limit, BMDL.
 *  external: Spec[]
 *  input:
 *   nparm is the number of parameters for Model B
 *   Anparm is the number of parameters for Model A
 *   xlk is the log-likelihood of the fitted model
 *   Dose is the BMD or upper dose limit
 *   pBak[] is the vector of fitted parameters for Model B
 *   pABak[] is the vector of fitted parameters for Model A
 *   D is a lower dose limit
 *   gtol is a small positive number (tolerance)
 *  output: lower confidence limit
 *****************************************************************/
//TODO: combine with BMDU_combofunc
double BMDL_combofunc(struct python_multitumor_analysis *pyAnal, struct python_multitumor_result *pyRes, double Dose, double D, double LR, double gtol, int *is_zero)
{ 	/* ck and LR are calculated in Multistage_ComboBMD() */

  int optite, nresm, *bind, CBnparm;
  int which, temprisk, lnParmMax;
  double fD, bmdl,  target, xmax, xlk2, xlk3, crisk;
  int i, j, nCall, k, nParmMax, ii = 0;
  double scale;

//  std::cout<<"inside BMDL_combofunc"<<std::endl;

  std::vector<double> adxmax;
  int nParms;
  fD = bmdl =  target = xmax = xlk2 = xlk3 = crisk = 0.0;
  optite = -5;
  nCall = 1;

  temprisk = 1;  //based on bmdparm.risk set to zero in original Multistage_ComboBMD

  adxmax.resize(pyRes->ndatasets, 0.0);



  CBnparm = lnParmMax = 0;
  for(i = 0; i < pyRes->ndatasets; i++)
  {
      int selModelIndex = pyRes->selectedModelIndex[i];
      struct python_dichotomous_analysis mod = pyAnal->models[i][selModelIndex];
      struct python_dichotomous_model_result modRes = pyRes->models[i][selModelIndex];
      xmax = mod.doses[0];
      CBnparm = CBnparm + modRes.nparms;
		
      if(modRes.nparms > lnParmMax){
         lnParmMax = modRes.nparms;
      }

//      std::cout<<"\n\nIn BMDL_combofunc, Tumor "<<i<<" data"<<std::endl;
//      std::cout<<"       DOSE     Inc    N"<<std::endl;
      for(j=0; j<mod.n; j++){
         if(mod.doses[j] > xmax){
            xmax = mod.doses[j];
         }
//         std::cout<<mod.doses[j]<<"   "<<mod.Y[j]<<"   "<<mod.n_group[j]<<std::endl;
      }
      adxmax[i] = xmax;
  }
  CBnparm = CBnparm + 1;

  /** rescale all doses to be: 0 <= Dose <= 1 **/
  xmax = adxmax[0];
  for(i=0;i<pyRes->ndatasets; i++)
    {
      if(adxmax[i] > xmax) xmax = adxmax[i];
    }
  scale = xmax;

  nParmMax = (int)lnParmMax;
  nParms = pyRes->ndatasets*nParmMax;
  std::vector<double> pdParms(nParms, 0.0);
  std::vector<double> pdParmsBak(nParms, 0.0);
  std::vector<double> pdParms2(nParms, 0.0);
  std::vector<double> pdVals(nParms, 0.0);
  std::vector<int> piSpec2(nParms, 0.0);


//  std::cout<<"\nIn BMDL_combofunc, pdParms[j](MLEs)"<<std::endl;
//  for(i=0; i<pyAnal->ndatasets; i++)
//    {
//      std::cout<<"Tumor "<<i<<"=>"<<std::endl;
//      for(j = 0; j<pyRes->models[i][0].nparms; j++)
//      {
//        std::cout<<pyRes->models[i][0].parms[j]<<"\t";
//      }
//      std::cout<<std::endl;
//    }

  k = -1;
  for(j=0; j<nParmMax; j++)
    {
      for(i=pyAnal->ndatasets-1; i>=0; i--)
	{
          k++;
          if(j<pyRes->models[i][0].nparms)
          {
            pdParms[k] = pyRes->models[i][0].parms[j];
            piSpec2[k] = 0.0;  //no user specified values
 	  }
	}
    }

//  std::cout<<"\n\nIn BMDL_combofunc, pdParms values (MLEs, k="<<k<<", nParms="<<nParms<<")"<<std::endl;
  i = 0;
//  for(j=0; j<pyAnal->ndatasets; j++){
//     std::cout<<"      Tumor "<<j<<"\t";
//  }
//  std::cout<<std::endl;

//  for(k = 0; k < nParmMax; k++)
//  {
//     for(j=0; j<pyAnal->ndatasets; j++){
//       std::cout<<pdParms[i++]<<"\t";
//     }
//     std::cout<<std::endl;
//   }

  j=0;
//  std::cout<<"\nIn BMDL_combofunc, Tumor Starting Values"<<std::endl;
//  for(i=0; i<pyAnal->ndatasets; i++)
//  {
//     std::cout<<"Tumor "<<i<<" => "<<pdParms[i]<<std::endl;
//  }

//  std::cout<<"\nMaximum Dose = "<<xmax<<std::endl;

  Dose = Dose/scale;
//  std::cout<<"Scale = "<<scale<<std::endl;

  which = 4;          /* Want a combined  lower confidence limit */

  target = (pyRes->combined_LL - LR);  /* The value we want the likelihood */
  
//  std::cout<< std::fixed << std::showpoint;
//  std::cout << std::setprecision(15);
//  std::cout<<"Combined Loglikelihood         "<<pyRes->combined_LL<<std::endl;
//  std::cout<<"Target                         "<<target<<std::endl;

  k = -1;

  for(j=0; j<nParmMax; j++)
    {
      for(i = pyRes->ndatasets-1; i>=0; i--)
      {
          int iParms = pyRes->models[i][0].nparms;

          k++;
          if (j <iParms){
            pdParmsBak[k] = pyRes->models[i][0].parms[j];
            pdParms[k] = pyRes->models[i][0].parms[j]*(pow(scale,(j)));
          } else {
            pdParmsBak[k] = pdParms[k] = BMDS_MISSING;
          }
      }
    }
//  /* One more step for the background terms */
  for(i=0; i<pyRes->ndatasets; i++){
      pdParms[i] = -log(1.0 - pdParms[i]);
  }
//  std::cout<<"\n\nValues BEFORE call "<<nCall<<" to getclmt_()";
//  std::cout<<"BMR="<<pyAnal->BMR<<" target="<<target<<std::endl;
//  std::cout<<"bmdl="<<bmdl<<" optite="<<optite<<std::endl;
  i = 0;
//  std::cout<<std::endl;
//  for(j=0; j<pyRes->ndatasets; j++){
//      std::cout<<"    Tumor "<<j<<"\t\t\t";
//  }
//  std::cout<<std::endl;
//  for(j=0; j<pyAnal->ndatasets; j++){
//      std::cout<<"Scaled | Unscaled\t\t";
//  }
//  std::cout<<std::endl;
//  for(k = 0; k < nParmMax; k++)
//    {
//      for(j=0; j<pyRes->ndatasets; j++){
//          std::cout<<pdParms[i]<<" | "<<pdParmsBak[i]<<"\t\t";
//	  i++;
//      }
//      std::cout<<std::endl;
//    }
//    std::cout<<"Dose(BMD)="<<Dose<<std::endl;
//    std::cout<<"log(BMD)="<<log(Dose)<<std::endl;

    double retVal;

    //bmdl calc
   

    bmdl = getclmt(pyAnal, pyRes, Dose, target, xmax, pdParms, true);
//    std::cout<<"getclmt returned bmdl = "<<bmdl<<std::endl;
//
//  fflush(fp_log);

//  fprintf(fp_log,"\n\nValues AFTER call %d to getclmt_()", nCall);
//  fprintf(fp_log,"BMR=%10.5g target=%10.5g\n",BMR, target);
//  fprintf(fp_log,"bmdl=%10.5g optite=%d", bmdl, optite);
//  i = 0;
//  fprintf(fp_log,"\n          ");
//  for(j = 1; j<=nT; j++)
//    fprintf(fp_log,"    Tumor %d\t", j);
//  fprintf(fp_log,"\n");
//
//  for(k = 0; k < nParmMax; k++)
//    {
//      for(j = 1; j <= nT; j++)
//	fprintf(fp_log,"%10.5g\t", pdParms[i++]);
//      fprintf(fp_log,"\n");
//    }
//  fflush(fp_log);
//  nCall++;
//
//  /* optite is a value that is passed back from GETCL which         */
//  /* determines whether the optimization was completed successfully */
//  /* If optite is less than 0, then it did not, and we want         */
//  /* to try a different starting point and recompute                */
//
//  if(optite < 0 )
//    {
//#ifdef MISC_OUT
//      /* Warn user */
//      fprintf(fp_out, "**** WARNING:  Completion code = %d.  Optimum not found. Trying new starting point****\n\n", optite);
//#endif
//      /* Try up to 10 times if needed */
//      for(ii = 0; ii < 10; ii++)
//	{
//#if !0
//	  GetNewParms2(pdParms, nParmMax);  /* Get a new starting point */
//#else
//	  /* Get original values */
//	  k = -1;
//	  for(i = 1; i <= nT; i++)
//	    {
//	      for(j = 1; j <= nParmMax; j++)
//		{
//		  k++;
//		  pdParms[k] = aParmList[i].pdParms[j] * (pow(scale,(j-1)));;
//		}
//	    }
//	  GetNewParms2(pdParms, nParmMax);  /* Get a new starting point */
//
//	  /* again, reparameterize p[0] */
//	  for(i = 0; i < nT; i++)
//	    {
//	      pdParms[i] = -log(1-aParmList[i+1].pdParms[1]);
//	    }
//
//#endif
//	  /* Try again */
//	  fprintf(fp_log,"\n\n\nValues BEFORE call %d to getclmt_()\n", nCall);
//	  fprintf(fp_log,"BMR=%10.5g target=%10.5g\n",BMR, target);
//	  fprintf(fp_log,"bmdl=%10.5g optite=%d", bmdl, optite);
//	  i = 0;
//	  fprintf(fp_log,"\n");
//	  for(j = 1; j<=nT; j++)
//	    fprintf(fp_log,"    Tumor %d\t", j);
//	  fprintf(fp_log,"\n");
//	  for(k = 0; k < nParmMax; k++)
//	    {
//	      for(j = 1; j <= nT; j++)
//		fprintf(fp_log,"%10.5g\t", pdParms[i++]);
//	      fprintf(fp_log,"\n");
//	    }
//	  fflush(fp_log);
//
//	  getclmt_(&which, &lnParmMax, &BMR, &Dose,
//		   &target, pdParms, piSpec2,
//		   pdParms, &temprisk, &bmdl,
//		   pdParms2, &optite, &nresm,
//		   bind, is_zero);
//
//	  fprintf(fp_log,"\n\nValues AFTER call %d to getclmt_()", nCall);
//	  fprintf(fp_log,"BMR=%10.5g target=%10.5g\n",BMR, target);
//	  fprintf(fp_log,"bmdl=%10.5g optite=%d", bmdl, optite);
//	  i = 0;
//	  fprintf(fp_log,"\n");
//	  for(j = 1; j<=nT; j++)
//	    fprintf(fp_log,"    Tumor %d\t", j);
//	  fprintf(fp_log,"\n");
//
//	  for(k = 0; k < nParmMax; k++)
//	    {
//	      for(j = 1; j <= nT; j++)
//		fprintf(fp_log,"%10.5g\t", pdParms[i++]);
//	      fprintf(fp_log,"\n");
//	    }
//	  nCall++;
//	  fflush(fp_log);
//
//	  /* if optite >= 0, it is successful, and we can stop */
//	  if(optite >= 0)
//	    break;
//#ifdef MISC_OUT
//	  /* otherwise, issues another warning, and continue trying */
//	  else
//	    fprintf(fp_out, "**** WARNING %d:  Completion code = %d trying new start****\n\n", ii, optite);
//#endif
//	} /* end for */
//
//    } /* end: if (optite < 0) */
//
//  if(optite < 0 )
//    {
//#ifdef MISC_OUT
//      /* Warn user */
//      fprintf(fp_out, "**** WARNING:  Completion code = %d.  Optimum not found. Trying new starting point****\n\n", optite);
//#endif
//      /* Try up to 10 times if needed */
//      for(ii = 0; ii < 10; ii++)
//	{
//	  /* Get original values */
//	  k = -1;
//	  for(i = 1; i <= nT; i++)
//	    {
//	      for(j = 1; j <= nParmMax; j++)
//		{
//		  k++;
//		  pdParms[k] = aParmList[i].pdParms[j] * (pow(scale,(j-1)));;
//		}
//	    }
//	  GetMoreParms2(pdParms, nParmMax);  /* Get a new starting point */
//
//	  /* again, reparameterize p[0] */
//	  for(i = 0; i < nT; i++)
//	    {
//	      pdParms[i] = -log(1-aParmList[i+1].pdParms[1]);
//	    }
//
//	  /* Try again */
//	  fprintf(fp_log,"\n\n\nValues BEFORE call %d to getclmt_()", nCall);
//	  fprintf(fp_log,"BMR=%10.5g target=%10.5g\n",BMR, target);
//	  fprintf(fp_log,"bmdl=%10.5g optite=%d", bmdl, optite);
//	  fprintf(fp_log,"\n");
//	  for(j = 1; j<=nT; j++)
//	    fprintf(fp_log,"    Tumor %d\t", j);
//	  fprintf(fp_log,"\n");
//	  i = 0;
//	  for(k = 0; k < nParmMax; k++)
//	    {
//	      for(j = 1; j <= nT; j++)
//		fprintf(fp_log,"%10.5g\t", pdParms[i++]);
//	      fprintf(fp_log,"\n");
//	    }
//	  fflush(fp_log);
//
//	  getclmt_(&which, &lnParmMax, &BMR, &Dose,
//		   &target, pdParms, piSpec2,
//		   pdParms, &temprisk, &bmdl,
//		   pdParms2, &optite, &nresm,
//		   bind, is_zero);
//
//	  fprintf(fp_log,"\n\nValues AFTER call %d to getclmt_()", nCall);
//	  fprintf(fp_log,"BMR=%10.5g target=%10.5g\n",BMR, target);
//	  fprintf(fp_log,"bmdl=%10.5g optite=%d", bmdl, optite);
//	  fprintf(fp_log,"\n");
//	  for(j = 1; j<=nT; j++)
//	    fprintf(fp_log,"    Tumor %d\t", j);
//	  fprintf(fp_log,"\n");
//	  i = 0;
//	  for(k = 0; k < nParmMax; k++)
//	    {
//	      for(j = 1; j <= nT; j++)
//		fprintf(fp_log,"%10.5g\t", pdParms[i++]);
//	      fprintf(fp_log,"\n");
//	    }
//	  nCall++;
//	  fflush(fp_log);
//
//	  /* if optite >= 0, it is successful, and we can stop */
//	  if(optite >= 0)
//	    break;
//#ifdef MISC_OUT
//	  /* otherwise, issues another warning, and continue trying */
//	  else
//	    fprintf(fp_out, "**** WARNING %d:  Completion code = %d trying new start****\n\n", ii, optite);
//#endif
//	} /* end for */
//
//    } /* end: if (optite < 0) */
//
//
//  /* Let user know if no optimum was found */
//  if(ii == 10)
//    {
//#ifdef MISC_OUT
//      fprintf(fp_out, "\nWarning:  completion code still negative");
//#endif
//      fprintf(fp_out, "\nBMDL did not converge for BMR = %f\n", BMR);
//      bmdl_bmr_flag = 1;
//    } /* end if */



//  std::cout<<"Here after bmdl calc"<<std::endl;
  pdParms2 = pdParms;
  int nT = 0;
  for(int i=0; i<pyRes->ndatasets; i++){
      nT++;
  }
  std::vector<std::vector<double>> ppdParms(nT, std::vector<double> (nParmMax, BMDS_MISSING));
//  std::cout<<"********** pdParms2 Values **********"<<std::endl;
//  for (int j=0; j<pyRes->ndatasets; j++){
//      std::cout<<"   Tumor " << j << "\t";
//  }
//  std::cout<<std::endl;

  k = -1;
  for (int j=0; j<nParmMax; j++)
  {
     for (int i=0; i<pyRes->ndatasets; i++)
     {
	  k++;
//          std::cout<<pdParms2[k]<<"\t";
          if (j < pyRes->models[i][0].nparms)
	  {
	    ppdParms[i][j] = pdParms2[k];
	    if(k < pyAnal->ndatasets)
            {                
	      ppdParms[i][j] = 1-exp(-pdParms2[k]);
	    }
	  }
     }
//     std::cout<<std::endl;
  }

  int flag = 1;
//  std::cout<<"scale="<<scale<<std::endl; 
  bmdl = exp(bmdl)*scale;
  
  xlk3 = ComboMaxLike2(flag,bmdl,&crisk, ppdParms, pyAnal, pyRes);
//  std::cout<<"BMDL_combofunc returns bmdl="<<bmdl<<std::endl;
  return bmdl;
}


double BMDU_combofunc(struct python_multitumor_analysis *pyAnal, struct python_multitumor_result *pyRes, double Dose, double D, double LR, double gtol, int *is_zero) {

  int optite, nresm, *bind, CBnparm;
  int which, temprisk, lnParmMax;
  double fD, bmdu,  target, xmax, xlk2, xlk3, crisk;
  int i, j, nCall, k, nParmMax, ii = 0;
  double scale;

//  std::cout<<"inside BMDU_combofunc"<<std::endl;

  std::vector<double> adxmax;
  int nParms;
  fD = bmdu =  target = xmax = xlk2 = xlk3 = crisk = 0.0;
  optite = -5;
  nCall = 1;

  temprisk = 1;  //based on bmdparm.risk set to zero in original Multistage_ComboBMD


  /* Get the degree of polynomial */
  adxmax.resize(pyRes->ndatasets, 0.0);

  CBnparm = lnParmMax = 0;
  for(i = 0; i < pyRes->ndatasets; i++)
  {
      int selModelIndex = pyRes->selectedModelIndex[i];
      struct python_dichotomous_analysis mod = pyAnal->models[i][selModelIndex];
      struct python_dichotomous_model_result modRes = pyRes->models[i][selModelIndex];
      xmax = mod.doses[0];
      CBnparm = CBnparm + modRes.nparms;
		
      if(modRes.nparms > lnParmMax){
         lnParmMax = modRes.nparms;
      }

//      std::cout<<"\n\nIn BMDU_combofunc, Tumor "<<i<<" data"<<std::endl;
//      std::cout<<"       DOSE     Inc    N"<<std::endl;
      for(j=0; j<mod.n; j++){
         if(mod.doses[j] > xmax){
            xmax = mod.doses[j];
         }
//         std::cout<<mod.doses[j]<<"   "<<mod.Y[j]<<"   "<<mod.n_group[j]<<std::endl;
      }
      adxmax[i] = xmax;
  }
  CBnparm = CBnparm + 1;

  /** rescale all doses to be: 0 <= Dose <= 1 **/
  xmax = adxmax[0];
  for(i=0;i<pyRes->ndatasets; i++)
    {
      if(adxmax[i] > xmax) xmax = adxmax[i];
    }
  scale = xmax;

  nParmMax = (int)lnParmMax;
  nParms = pyRes->ndatasets*nParmMax;
  std::vector<double> pdParms(nParms, 0.0);
  std::vector<double> pdParmsBak(nParms, 0.0);
  std::vector<double> pdParms2(nParms, 0.0);
  std::vector<double> pdVals(nParms, 0.0);
  std::vector<int> piSpec2(nParms, 0.0);


//  std::cout<<"\nIn BMDU_combofunc, pdParms[j](MLEs)"<<std::endl;
//  for(i=0; i<pyAnal->ndatasets; i++)
//    {
//      std::cout<<"Tumor "<<i<<"=>"<<std::endl;
//      for(j = 0; j<pyRes->models[i][0].nparms; j++)
//      {
//        std::cout<<pyRes->models[i][0].parms[j]<<"\t";
//      }
//      std::cout<<std::endl;
//    }

  k = -1;
  for(j=0; j<nParmMax; j++)
    {
      for(i=pyAnal->ndatasets-1; i>=0; i--)
	{
          k++;
          if(j<pyRes->models[i][0].nparms)
          {
            pdParms[k] = pyRes->models[i][0].parms[j];
            piSpec2[k] = 0.0;  //no user specified values
 	  }
	}
    }

//  std::cout<<"\n\nIn BMDU_combofunc, pdParms values (MLEs, k="<<k<<", nParms="<<nParms<<")"<<std::endl;
  i = 0;
//  for(j=0; j<pyAnal->ndatasets; j++){
//     std::cout<<"      Tumor "<<j<<"\t";
//  }
//  std::cout<<std::endl;

//  for(k = 0; k < nParmMax; k++)
//  {
//     for(j=0; j<pyAnal->ndatasets; j++){
//       std::cout<<pdParms[i++]<<"\t";
//     }
//     std::cout<<std::endl;
//   }

  j=0;
//  std::cout<<"\nIn BMDU_combofunc, Tumor Starting Values"<<std::endl;
//  for(i=0; i<pyAnal->ndatasets; i++)
//  {
//     std::cout<<"Tumor "<<i<<" => "<<pdParms[i]<<std::endl;
//  }
//
//  std::cout<<"\nMaximum Dose = "<<xmax<<std::endl;

  Dose = Dose/scale;
//  std::cout<<"Scale = "<<scale<<std::endl;

  which = 5;          /* Want a combined  upper confidence limit */

  target = (pyRes->combined_LL - LR);  /* The value we want the likelihood */
//  std::cout<<"LR = "<<LR<<std::endl;
//  std::cout<<"Combined Loglikelihood         "<<pyRes->combined_LL<<std::endl;
//  std::cout<<"Target                         "<<target<<std::endl;

  k = -1;

  for(j=0; j<nParmMax; j++)
    {
      for(i = pyRes->ndatasets-1; i>=0; i--)
      {
          int iParms = pyRes->models[i][0].nparms;

          k++;
          if (j <iParms){
            pdParmsBak[k] = pyRes->models[i][0].parms[j];
            pdParms[k] = pyRes->models[i][0].parms[j]*(pow(scale,(j)));
          } else {
            pdParmsBak[k] = pdParms[k] = BMDS_MISSING;
          }
      }
    }
  for(i=0; i<pyRes->ndatasets; i++){
      pdParms[i] = -log(1.0 - pdParms[i]);
  }
//  std::cout<<"\n\nValues BEFORE call "<<nCall<<" to getclmt_()";
//  std::cout<<"BMR="<<pyAnal->BMR<<" target="<<target<<std::endl;
//  std::cout<<"bmdu="<<bmdu<<" optite="<<optite<<std::endl;
  i = 0;
//  std::cout<<std::endl;
//  for(j=0; j<pyRes->ndatasets; j++){
//      std::cout<<"    Tumor "<<j<<"\t\t\t";
//  }
//  std::cout<<std::endl;
//  for(j=0; j<pyAnal->ndatasets; j++){
//      std::cout<<"Scaled | Unscaled\t\t";
//  }
//  std::cout<<std::endl;
//  for(k = 0; k < nParmMax; k++)
//    {
//      for(j=0; j<pyRes->ndatasets; j++){
//          std::cout<<pdParms[i]<<" | "<<pdParmsBak[i]<<"\t\t";
//	  i++;
//      }
//      std::cout<<std::endl;
//    }
//    std::cout<<"Dose(BMD)="<<Dose<<std::endl;
//    std::cout<<"log(BMD)="<<log(Dose)<<std::endl;

    double retVal;

    bmdu = getclmt(pyAnal, pyRes, Dose, target, xmax, pdParms, false);
//  fflush(fp_log);

//  getclmt_(&which, &lnParmMax, &BMR, &Dose,
//	   &target, pdParms, piSpec2,
//	   pdParms, &temprisk, &bmdl,
//	   pdParms2, &optite, &nresm,
//	   bind, is_zero);
//
//  fprintf(fp_log,"\n\nValues AFTER call %d to getclmt_()", nCall);
//  fprintf(fp_log,"BMR=%10.5g target=%10.5g\n",BMR, target);
//  fprintf(fp_log,"bmdl=%10.5g optite=%d", bmdl, optite);
//  i = 0;
//  fprintf(fp_log,"\n          ");
//  for(j = 1; j<=nT; j++)
//    fprintf(fp_log,"    Tumor %d\t", j);
//  fprintf(fp_log,"\n");
//
//  for(k = 0; k < nParmMax; k++)
//    {
//      for(j = 1; j <= nT; j++)
//	fprintf(fp_log,"%10.5g\t", pdParms[i++]);
//      fprintf(fp_log,"\n");
//    }
//  fflush(fp_log);
//  nCall++;
//
//  /* optite is a value that is passed back from GETCL which         */
//  /* determines whether the optimization was completed successfully */
//  /* If optite is less than 0, then it did not, and we want         */
//  /* to try a different starting point and recompute                */
//
//  if(optite < 0 )
//    {
//#ifdef MISC_OUT
//      /* Warn user */
//      fprintf(fp_out, "**** WARNING:  Completion code = %d.  Optimum not found. Trying new starting point****\n\n", optite);
//#endif
//      /* Try up to 10 times if needed */
//      for(ii = 0; ii < 10; ii++)
//	{
//#if !0
//	  GetNewParms2(pdParms, nParmMax);  /* Get a new starting point */
//#else
//	  /* Get original values */
//	  k = -1;
//	  for(i = 1; i <= nT; i++)
//	    {
//	      for(j = 1; j <= nParmMax; j++)
//		{
//		  k++;
//		  pdParms[k] = aParmList[i].pdParms[j] * (pow(scale,(j-1)));;
//		}
//	    }
//	  GetNewParms2(pdParms, nParmMax);  /* Get a new starting point */
//
//	  /* again, reparameterize p[0] */
//	  for(i = 0; i < nT; i++)
//	    {
//	      pdParms[i] = -log(1-aParmList[i+1].pdParms[1]);
//	    }
//
//#endif
//	  /* Try again */
//	  fprintf(fp_log,"\n\n\nValues BEFORE call %d to getclmt_()\n", nCall);
//	  fprintf(fp_log,"BMR=%10.5g target=%10.5g\n",BMR, target);
//	  fprintf(fp_log,"bmdl=%10.5g optite=%d", bmdl, optite);
//	  i = 0;
//	  fprintf(fp_log,"\n");
//	  for(j = 1; j<=nT; j++)
//	    fprintf(fp_log,"    Tumor %d\t", j);
//	  fprintf(fp_log,"\n");
//	  for(k = 0; k < nParmMax; k++)
//	    {
//	      for(j = 1; j <= nT; j++)
//		fprintf(fp_log,"%10.5g\t", pdParms[i++]);
//	      fprintf(fp_log,"\n");
//	    }
//	  fflush(fp_log);
//
//	  getclmt_(&which, &lnParmMax, &BMR, &Dose,
//		   &target, pdParms, piSpec2,
//		   pdParms, &temprisk, &bmdl,
//		   pdParms2, &optite, &nresm,
//		   bind, is_zero);
//
//	  fprintf(fp_log,"\n\nValues AFTER call %d to getclmt_()", nCall);
//	  fprintf(fp_log,"BMR=%10.5g target=%10.5g\n",BMR, target);
//	  fprintf(fp_log,"bmdl=%10.5g optite=%d", bmdl, optite);
//	  i = 0;
//	  fprintf(fp_log,"\n");
//	  for(j = 1; j<=nT; j++)
//	    fprintf(fp_log,"    Tumor %d\t", j);
//	  fprintf(fp_log,"\n");
//
//	  for(k = 0; k < nParmMax; k++)
//	    {
//	      for(j = 1; j <= nT; j++)
//		fprintf(fp_log,"%10.5g\t", pdParms[i++]);
//	      fprintf(fp_log,"\n");
//	    }
//	  nCall++;
//	  fflush(fp_log);
//
//	  /* if optite >= 0, it is successful, and we can stop */
//	  if(optite >= 0)
//	    break;
//#ifdef MISC_OUT
//	  /* otherwise, issues another warning, and continue trying */
//	  else
//	    fprintf(fp_out, "**** WARNING %d:  Completion code = %d trying new start****\n\n", ii, optite);
//#endif
//	} /* end for */
//
//    } /* end: if (optite < 0) */
//
//  if(optite < 0 )
//    {
//#ifdef MISC_OUT
//      /* Warn user */
//      fprintf(fp_out, "**** WARNING:  Completion code = %d.  Optimum not found. Trying new starting point****\n\n", optite);
//#endif
//      /* Try up to 10 times if needed */
//      for(ii = 0; ii < 10; ii++)
//	{
//	  /* Get original values */
//	  k = -1;
//	  for(i = 1; i <= nT; i++)
//	    {
//	      for(j = 1; j <= nParmMax; j++)
//		{
//		  k++;
//		  pdParms[k] = aParmList[i].pdParms[j] * (pow(scale,(j-1)));;
//		}
//	    }
//	  GetMoreParms2(pdParms, nParmMax);  /* Get a new starting point */
//
//	  /* again, reparameterize p[0] */
//	  for(i = 0; i < nT; i++)
//	    {
//	      pdParms[i] = -log(1-aParmList[i+1].pdParms[1]);
//	    }
//
//	  /* Try again */
//	  fprintf(fp_log,"\n\n\nValues BEFORE call %d to getclmt_()", nCall);
//	  fprintf(fp_log,"BMR=%10.5g target=%10.5g\n",BMR, target);
//	  fprintf(fp_log,"bmdl=%10.5g optite=%d", bmdl, optite);
//	  fprintf(fp_log,"\n");
//	  for(j = 1; j<=nT; j++)
//	    fprintf(fp_log,"    Tumor %d\t", j);
//	  fprintf(fp_log,"\n");
//	  i = 0;
//	  for(k = 0; k < nParmMax; k++)
//	    {
//	      for(j = 1; j <= nT; j++)
//		fprintf(fp_log,"%10.5g\t", pdParms[i++]);
//	      fprintf(fp_log,"\n");
//	    }
//	  fflush(fp_log);
//
//	  getclmt_(&which, &lnParmMax, &BMR, &Dose,
//		   &target, pdParms, piSpec2,
//		   pdParms, &temprisk, &bmdl,
//		   pdParms2, &optite, &nresm,
//		   bind, is_zero);
//
//	  fprintf(fp_log,"\n\nValues AFTER call %d to getclmt_()", nCall);
//	  fprintf(fp_log,"BMR=%10.5g target=%10.5g\n",BMR, target);
//	  fprintf(fp_log,"bmdl=%10.5g optite=%d", bmdl, optite);
//	  fprintf(fp_log,"\n");
//	  for(j = 1; j<=nT; j++)
//	    fprintf(fp_log,"    Tumor %d\t", j);
//	  fprintf(fp_log,"\n");
//	  i = 0;
//	  for(k = 0; k < nParmMax; k++)
//	    {
//	      for(j = 1; j <= nT; j++)
//		fprintf(fp_log,"%10.5g\t", pdParms[i++]);
//	      fprintf(fp_log,"\n");
//	    }
//	  nCall++;
//	  fflush(fp_log);
//
//	  /* if optite >= 0, it is successful, and we can stop */
//	  if(optite >= 0)
//	    break;
//#ifdef MISC_OUT
//	  /* otherwise, issues another warning, and continue trying */
//	  else
//	    fprintf(fp_out, "**** WARNING %d:  Completion code = %d trying new start****\n\n", ii, optite);
//#endif
//	} /* end for */
//
//    } /* end: if (optite < 0) */
//
//
//  /* Let user know if no optimum was found */
//  if(ii == 10)
//    {
//#ifdef MISC_OUT
//      fprintf(fp_out, "\nWarning:  completion code still negative");
//#endif
//      fprintf(fp_out, "\nBMDL did not converge for BMR = %f\n", BMR);
//      bmdl_bmr_flag = 1;
//    } /* end if */


//  double **ppdParms;
//  ppdParms = DMATRIX (1, nT, 1, nParmMax);

  pdParms2 = pdParms;
  int nT = 0;
  for(int i=0; i<pyRes->ndatasets; i++){
      nT++;
  }
  std::vector<std::vector<double>> ppdParms(nT, std::vector<double> (nParmMax, BMDS_MISSING));
//  std::cout<<"********** pdParms2 Values **********"<<std::endl;
//  for (int j=0; j<pyRes->ndatasets; j++){
//      std::cout<<"   Tumor " << j << "\t";
//  }
//  std::cout<<std::endl;

  k = -1;
  for (int j=0; j<nParmMax; j++)
  {
     for (int i=0; i<pyRes->ndatasets; i++)
     {
	  k++;
//          std::cout<<pdParms2[k]<<"\t";
          if (j < pyRes->models[i][0].nparms)
	  {
	    ppdParms[i][j] = pdParms2[k];
	    if(k < pyAnal->ndatasets)
            {                
	      ppdParms[i][j] = 1-exp(-pdParms2[k]);
	    }
	  }
     }
//     std::cout<<std::endl;
  }

  int flag = 1; 
  bmdu = exp(bmdu)*scale;
  
  xlk3 = ComboMaxLike2(flag,bmdu,&crisk, ppdParms, pyAnal, pyRes);
  return bmdu;
}


double ComboMaxLike2(int flag, double dose, double *crisk, std::vector<std::vector<double>> p, python_multitumor_analysis *pyAnal, python_multitumor_result *pyRes){

  int nObs, nT, nParms;
  double prob, like, dSumParm1, pr, bkg;

  prob = like = dSumParm1 = 0.0;

  for(int n=0; n<pyRes->ndatasets; n++){
      dSumParm1 += p[n][0];
      nObs = pyAnal->models[n][0].n;
      int selIndex = pyRes->selectedModelIndex[n];
      nParms = pyRes->models[n][selIndex].nparms;
      for (int i=0; i<nObs; i++){
         double D, Yn, Yp;
         D = pyAnal->models[n][selIndex].doses[i];
         Yp = pyAnal->models[n][selIndex].Y[i];
         Yn = pyAnal->models[n][selIndex].n_group[i] - Yp;
         prob = p[n][nParms-1];
         if(n==0 && flag==1 && nParms==2){
           prob = p[0][1];
           for(int nt=1; nt<pyRes->ndatasets; nt++){
               prob -= p[nt][1];
           }
         }

         for (int j=nParms-1; j>=0; j--){
           if (n==0 && flag==1 && j==1){
              pr = p[0][1];
              for (int nt=1; nt<pyRes->ndatasets; nt++){
                    pr -= p[nt][1];
              }
              prob = D*prob + pr;
           } else {
              prob = D*prob + p[n][j]; 
           }
         }
         prob = (1-exp(-1.0* prob));
         if ((prob==0) || (prob == 1)){
           if(Yp <=0 || Yn <=0){
             like += 0;
           } else {
             if (prob == 1){
               like += Yn*(-710);
             } else {
               like += Yp*(-710);
             }
           }
         } else {
           like += Yp*log(prob) + Yn*log(1-prob);
         }
      }
  }


  bkg = 1.0 - exp(-1.0*(dSumParm1));


  for(int n=0; n<pyRes->ndatasets; n++){
    int selIndex = pyRes->selectedModelIndex[n];
    nParms = pyRes->models[n][selIndex].nparms;

    prob = p[n][nParms-1];
    if (n ==0 && flag == 1 && nParms == 2){
       prob = p[0][1];
       for(int nt = 1; nt<pyRes->ndatasets; nt++){
           prob -= p[nt][1];
       }
    }
    for (int j=nParms-1; j>=0; j--){
       if (n==0 && flag == 1 && j==1){
         pr = p[0][1];
         for (int nt = 1; nt<pyRes->ndatasets; nt++){
             pr -= p[nt][1];
         }
         prob = dose*prob + pr;
       } else {
         prob = dose*prob + p[n][j];
       }
    }

  }

  if (bkg == 1.0){
    *crisk = 0.0;
  } else {
    *crisk = ((1.0 - exp(-1.0*(prob))) - bkg)/(1.0-bkg);
  }

  return like;
}


/***************************************************************************
 * Multistage_ComboBMD -- Used to calculate the BMD and BMDL for combined 
 *                         Multistage models A and B
 *  external: bmdparm
 *  input:      *   nparm is the number of parameters
 *   nparm is the number of parameters/
 *   p[] is the vector of fitted parameters
 *   gtol is a small positive number
 *   iter is not used ????????
 *   xlk is the sum of log-likelihood for the fitted models (A + B)
 *   Rlevel[] is the vector of BMR's
 *   Bmdl[] is the vector of BMDL's for the BMR's
 *   Bmdu[] is the vector of BMDU's for the BMR's
 *   BMD is the benchmark dose
 *  output: BMD, Bmdl[], prints BMDL
 ****************************************************************************/
void Multistage_ComboBMD (struct python_multitumor_analysis *pyAnal, struct python_multitumor_result *pyRes){

   int cnparm, selIndex, nT, is_zero;
   double LR, xa, xb, D, fa, fb, Drange, cxmax, ck, poly;
   double gtol = 1e-12;
   double tol;  //for zeroin function
   std::vector<double> cp;


   nT = pyRes->ndatasets;
   //find largest nparm and largest dose
   cnparm = 0;
   cxmax = 0;
   for(int i=0; i<nT; i++){
        selIndex = pyRes->selectedModelIndex[i];
        if(pyRes->models[i][selIndex].nparms > cnparm){
           cnparm = pyRes->models[i][selIndex].nparms;
        }
        double tmpMax = *max_element(std::begin(pyAnal->models[i][selIndex].doses), std::end(pyAnal->models[i][selIndex].doses));
        if(cxmax < tmpMax) cxmax = tmpMax;
   }
   //add all model parameters to combined p
   cp.resize(cnparm, 0.0);
   for(int i=0; i<nT; i++) {
       selIndex = pyRes->selectedModelIndex[i];
       cp[0] = cp[0] - log(1.0 - pyRes->models[i][selIndex].parms[0]);
       for(int j=1; j<pyRes->models[i][selIndex].nparms; j++){
          cp[j] = cp[j] + pyRes->models[i][selIndex].parms[j];  
       }
   }
   //compute chi-squared value
   double cl = 1.0 - pyAnal->alpha;
   if (cl<0.5){
     LR = 0.5*gsl_cdf_chisq_Pinv(1.0-2.0*cl, 1);
   } else {
     LR = 0.5*gsl_cdf_chisq_Pinv(2.0*cl-1.0, 1);
   }

   ck = -log(1-pyAnal->BMR);  //Extra risk
   //solve the BMD
   xa = D = 0.0;
   fa = -ck;  //note:  ck>0.0
   fb = fa;
   Drange = cxmax; 

   int k=1;
   while(k<300 && fb<0){
      fa=fb;
      xa=D;
      D=Drange*k/100.0;
      poly=cp[cnparm-1];
      for (int j=cnparm-2; j>0; j--){
         poly=poly*D+cp[j];
      }
      poly=poly*D;
      fb=poly-ck;
      k++;
   }

   if (fb<0) std::cout<<"BMD Computation failed.  BMD is larger than three times maximum input doses."<<std::endl;
   xb = D;
   tol = 1.0e-15;

   //compute BMD
   //BMD_func works on log scale, so convert xa and xb to logs
   if (xa==0.0) xa = -690.7755;
   else xa = log(xa);
   xb = log(xb);

   xb = zeroin(xa, xb, tol, BMD_func, cnparm, &cp[0], ck);
   xa = exp(xa);
   xb = exp(xb);
   pyRes->BMD = xb;

   is_zero = 0;

   pyRes->BMDL = BMDL_combofunc(pyAnal, pyRes, xb, xa, LR, tol, &is_zero);
//   std::cout<<"pyRes->BMDL="<<pyRes->BMDL<<std::endl;

   pyRes->BMDU = BMDU_combofunc(pyAnal, pyRes, xb, xa, LR, tol, &is_zero);
//   std::cout<<"pyRes->BMDU="<<pyRes->BMDU<<std::endl;

}







double objfunc_bmdl(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data){
   //obj function of form F(BMD, beta) = BMD, where BMD=X[0] and beta=X[i] where i>0
   if (!grad.empty()){
      //fill all gradients to zero (grad[1] to grad[n] are all zero
      // because objective function only depends on X[0] 
      std::fill(grad.begin(), grad.end(), 0);
      //set first grad to 1, since x[1] should be the BMD (Dose)
      grad[0] = 1.0;
   }
   return x[0];  
}


double objfunc_bmdu(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data){
   //obj function of form F(BMD, beta) = BMD, where BMD=X[0] and beta=X[i] where i>0
   if (!grad.empty()){
      //fill all gradients to zero (grad[1] to grad[n] are all zero
      // because objective function only depends on X[0] 
      std::fill(grad.begin(), grad.end(), 0);
      //set first grad to 1, since x[1] should be the BMD (Dose)
      grad[0] = -1.0;
   }
   return -1*x[0];  
}

double myEqualityConstraint(const std::vector<double> &x, std::vector<double> &grad, void *data){
  
   msComboEq *d = reinterpret_cast<msComboEq*>(data);
   double bmr = d->bmr;
   int nT = d->nT;
   std::vector<int> degree = d->degree;

   double D = exp(x[0]);
   int iIndex = x.size() - 1;
   double sum2 = 0.0;
   double sum = 0.0;

   if (!grad.empty()){
     for (int l=nT-1; l>=0; l--){
        sum = 0.0;
        for (int k=degree[l]; k>0; k--){
           sum = sum*D + k*x[iIndex]*D;
           iIndex -= 1;
        }
        iIndex -= 1;
        sum2 += sum;
     } 
     grad[0] = sum2;

     iIndex = 1;
     for (int k=0; k<nT; k++){
	for (int j=0; j<=degree[k]; j++){
	   if (j==0){
	      grad[iIndex]=0;
           } else {
	      grad[iIndex] = pow(D,j);
	   }
	   iIndex +=1;
	}

     }
   }

   //equality constraint calc
   sum = log(1.0 - bmr);
   iIndex = x.size() - 1;
   double sum3 = 0.0;

   for (int l=nT-1; l>=0; l--){
     sum2 = 0.0;
     for (int k=degree[l]; k>0; k--){
        sum2 = sum2*D + x[iIndex]*D;
        iIndex -= 1;    
     }
     iIndex -= 1;
     sum3 += sum2;
   }

   return sum + sum3;
}


double myInequalityConstraint1(const std::vector<double> &x, std::vector<double> &grad, void *data){

   msComboInEq *d = reinterpret_cast<msComboInEq*>(data);
   double target = d->target;
   int nT = d->nT;
   std::vector<int> nObs = d->nObs;
   std::vector<int> degree = d->degree;
   std::vector<std::vector<double>> doses = d->doses;
   std::vector<std::vector<double>> Y = d->Y;
   std::vector<std::vector<double>> n_group = d->n_group;

   int m = nT;
   int iOffset = 1;
   double resid = 0.0; 
   if (!grad.empty()){
     for (size_t i=0; i<grad.size(); i++){
        grad[i] = 0.0;
     }
     for (int l=0; l<nT; l++){
       m = m-1;
       double iTop = degree[l] + iOffset;
       double iBottom = iOffset;
       for (int k=0; k<nObs[l]; k++){
         double sum = x[iTop];
         for (int j=iTop-1; j>=iBottom; j--){
            sum = sum * doses[m][k] + x[j];
         }
         if (sum < 0) sum = 0.0;
         double P = 1.0 - exp(-1.0*sum);
         resid = (Y[m][k]*dslog(P) - (n_group[m][k]-Y[m][k])*dslog(1.0-P))*(P-1.0);
         int iIndex = iTop-1;
         for(int j=degree[l]-1; j>=0; j--){
           grad[iIndex] = grad[iIndex] + resid*(pow(doses[m][k],j));
           iIndex = iIndex-1;
         }
       }  
       iOffset = iTop + 1; 
     }
   }


   double sum2 = 0;
   iOffset = 1;
   m = nT;
   for (int l=0; l<nT; l++){
     m = m-1;
     double iTop = degree[l] + iOffset;
     double iBottom = iOffset;
     for (int k=0; k<nObs[l]; k++){
       double sum = x[iTop];
       for (int j= iTop-1; j>=iBottom; j--){
          sum = sum*doses[m][k] + x[j];
       }
       if (sum < 0) sum = 0.0;
       double P = 1.0 - exp(-1*sum);
       sum2 = sum2 + Y[m][k]*slog(P) + (n_group[m][k]-Y[m][k]) * slog(1.0 - P);
     }
     iOffset = iTop + 1;
   }
   return  target - sum2;

}


double slog(double X){
   double coefs[4] = {6.7165863851209542e50,-2.0154759155362862e+35, 2.0169759155362859e+19,-710};
   if (X >= 1e-16){
     return log(X);
   } else {
     double v = 0.0;
     for (int i=0; i<4; i++){
       v = X * v + coefs[i];
     }
     return v;
   }
}

double dslog(double P){
   if (P >= 1e-10){
     return 1.0/P;
   } else {
     return 2.0e10 - 1.0e20 * P;
   }
}



void BMDS_ENTRY_API __stdcall runBMDSDichoAnalysis(struct dichotomous_analysis *anal, struct dichotomous_model_result *res, struct dichotomous_GOF *gof, struct BMDS_results *bmdsRes, struct dicho_AOD *bmdsAOD){
 
//  std::cout<<"degree = "<<anal->degree;
//  std::cout<<"degree:"<<anal->degree<<std::endl;
//        for (int k=0; k<anal->prior_cols*anal->parms; k++){
//          std::cout<<anal->prior[k]<<", ";
//        }
//        std::cout<<std::endl;
 
  bmdsRes->validResult = false;
  bmdsRes->slopeFactor = BMDS_MISSING;
  bmdsRes->BMD = BMDS_MISSING;
  bmdsRes->BMDL = BMDS_MISSING;
  bmdsRes->BMDU = BMDS_MISSING;
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


  bmdsRes->BMD_MA = BMDS_MISSING;
  bmdsRes->BMDL_MA = BMDS_MISSING;
  bmdsRes->BMDU_MA = BMDS_MISSING;

  estimate_ma_laplace_dicho(MA, DA, res);


  collect_dichoMA_bmd_values(MA, res, bmdsRes, DA->alpha);
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

  bmdsRes->BMD = BMDS_MISSING;
  bmdsRes->BMDL = BMDS_MISSING;
  bmdsRes->BMDU = BMDS_MISSING;
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
//    gof->n = GOFanal.n;
  }
  gof->n = GOFanal.n;
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
  anal->isIncreasing = pyAnal->isIncreasing;
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

  bool validated = false;
  if (pyAnal->suff_stat){
    validated = pyAnal->n == pyAnal->doses.size() && pyAnal->doses.size() == pyAnal->Y.size() && pyAnal->doses.size() == pyAnal->n_group.size();
  } else {
    validated = pyAnal->n == pyAnal->doses.size() && pyAnal->doses.size() == pyAnal->Y.size();
  }

  if (validated){
    for (int i=0; i<pyAnal->n; i++){
      anal->Y[i] = pyAnal->Y[i];
      anal->doses[i] = pyAnal->doses[i];
    }
  }
  if (validated && pyAnal->suff_stat){
    for (int i=0; i<pyAnal->n; i++){
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
     //only all individual multistage models if degree == 0
     //else only run the specified model
     if (pyAnal->degree[i] ==0){
       for (int j=0; j<pyAnal->nmodels[i]; j++){
         pythonBMDSDicho(&pyAnal->models[i][j], &pyRes->models[i][j]);  
       } 
     } else {
       pythonBMDSDicho(&pyAnal->models[i][0], &pyRes->models[i][0]);   
     }
   }

   //select models
   selectMultitumorModel(pyAnal, pyRes);

   //create new pyAnal and pyRes only containing selected models
   //Note: ndatasets may differ from previous structs (pyAnal & pyRes) depending on whether any datasets were rejected
   struct python_multitumor_analysis anal;
   struct python_multitumor_result res;


   anal.n = pyAnal->n;
   anal.BMR = pyAnal->BMR;
   anal.alpha = pyAnal->alpha;
   anal.prior_cols = pyAnal->prior_cols;
   pyRes->validResult.clear();
   for (int dataset=0; dataset<pyAnal->ndatasets; dataset++){
     int selIndex = pyRes->selectedModelIndex[dataset];
     if (selIndex >= 0){
        pyRes->validResult.push_back(true);
        std::vector<python_dichotomous_analysis> modGroup;
        std::vector<python_dichotomous_model_result> modResGroup;
        struct python_dichotomous_analysis modAnal;
        struct python_dichotomous_model_result modRes;
        
        anal.nmodels.push_back(1);
        modAnal = pyAnal->models[dataset][selIndex];
        modGroup.push_back(modAnal);
	anal.degree.push_back(modAnal.degree);
	anal.models.push_back(modGroup);
	anal.n.push_back(modAnal.n);
	anal.ndatasets++;

	res.nmodels.push_back(1);
	res.selectedModelIndex.push_back(0);  //selected index is always zero since only one model
        modRes = pyRes->models[dataset][selIndex];
        modResGroup.push_back(modRes);
	res.models.push_back(modResGroup);

     } else {
	//std::cout<<"Warning: multistage model not selected for dataset:"<<dataset<<std::endl;
	//std::cout<<"This dataset will be dropped."<<std::endl;
	pyRes->validResult.push_back(false);
     }
   }
   anal.ndatasets = anal.nmodels.size();
   res.ndatasets = anal.ndatasets;

//   std::cout<<"running multitumor model with "<<anal.ndatasets<<" datasets"<<std::endl;
   //run MSCombo
   runMultitumorModel(&anal, &res);  

   pyRes->BMD = res.BMD;
   pyRes->BMDL = res.BMDL;
   pyRes->BMDU = res.BMDU;
   pyRes->slopeFactor = res.slopeFactor;
   pyRes->combined_LL = res.combined_LL;
   pyRes->combined_LL_const = res.combined_LL_const;
}  


void selectMultitumorModel(struct python_multitumor_analysis *pyAnal, struct python_multitumor_result *pyRes){


   if (pyRes->selectedModelIndex.size() > 0){
      pyRes->selectedModelIndex.clear();
   }

   for (int i=0; i<pyAnal->ndatasets; i++){ 
      if(pyAnal->degree[i] == 0){
        int selectedIndex = selectBestMultitumorModel(pyAnal->models[i], pyRes->models[i]);
        pyRes->selectedModelIndex.push_back(selectedIndex);
      } else {
	//in this case the only the user selected model was run
	pyRes->selectedModelIndex.push_back(0);
      }
   }

}

int selectBestMultitumorModel(std::vector<python_dichotomous_analysis> &analModels, std::vector<python_dichotomous_model_result> &resModels){

   double bmdDoseSRLimit = 2.0; // scaled residual at BMD dose < value
   double controlDoseSRLimit = 2.0; // scaled residual at control dose < value
   double pValueMin = 0.05; // p-value > value

   std::vector<bool> modelHitBoundary;  //checks whether specific model hit boundary
   std::vector<bool> adequateFit;
   bool anyHitBoundary = false;
   bool anyAdequateFit = false;

   for (int i=0; i<analModels.size(); i++){
      //check parameter boundaries
      modelHitBoundary.push_back(false);
      std::vector<bool> bounded = resModels[i].bmdsRes.bounded;
      for (int j=0; j<bounded.size(); j++){
         if (bounded[j]){
            modelHitBoundary[i] = true;
//	    std::cout<<"    hit boundary for parm:"<<j<<std::endl;
	 }
      }
      //check adequate fit
      //1. scaled residual at BMD dose < bmdDoseSRLimit
      //2. scaled residual at control dose < controlDoseSRLimit
      //3. If 1 & 2 are met, p-value > pValueMin
      //Note: these values are all currently hardcoded.  Need to account for/allow user changes
      adequateFit.push_back(false);
      if (resModels[i].getSRAtDose(resModels[i].bmdsRes.BMD, analModels[i].doses) < bmdDoseSRLimit && resModels[i].gof.residual[0] < controlDoseSRLimit && resModels[i].gof.p_value > pValueMin){
        adequateFit[i] = true;
      }
   }


   //first only examine models up to k-2
   anyHitBoundary = false;
   anyAdequateFit = false;
   for (int i=0; i<analModels.size()-1; i++){
     if(modelHitBoundary[i]){
        anyHitBoundary = true;
     }
     if(adequateFit[i]){
        anyAdequateFit = true;
     }
   }

   if (anyHitBoundary){
//      std::cout<<"Using SOP step 2 due to model hitting boundary"<<std::endl;
      //at least one model had a parameter that hit a boundary
      //step 2 under SOP instructions
      if(adequateFit[0] && adequateFit[1]){
         //both linear and quadratic models fit
         if(!modelHitBoundary[0] && !modelHitBoundary[1]){
            //choose model with lowest AIC
            if (round_to(resModels[1].bmdsRes.AIC, BMDS_EPS) < round_to(resModels[0].bmdsRes.AIC, BMDS_EPS)){
	       return 1;
	    } else {
	       return 0;
	    }
	 } else {
           //choose model with lowest BMDL
	   if (round_to(resModels[1].bmdsRes.BMDL, BMDS_EPS) < round_to(resModels[0].bmdsRes.BMDL, BMDS_EPS)){
              return 1;
	   } else {
	      return 0;
	   }
	 }	 
      } else if(adequateFit[0]){
        //only linear fits
	return 0;
      } else if(adequateFit[1]){
        //only quadratic fits
	return 1;
      } else {
	//neither fit - refer to SWG
	return -1;
      }
      //return 1;

   } else {
      //no models had a paramter hit the boundary
      //step 1 under SOP instructions
      if(anyAdequateFit){
         //use lowest AIC model that adequately fits up to k-2 model (exclude k-1 model)
         std::vector<double> modelAIC;
	 for (int i=0; i<resModels.size()-1; i++){
            modelAIC.push_back(resModels[i].bmdsRes.AIC);  
	 }
	 //now find the lowest AIC
         int index = std::distance(std::begin(modelAIC), std::min_element(std::begin(modelAIC), std::end(modelAIC)));
         return index;
      } else {
	 //no adequateFit models.  Examine k-1 model
         if(adequateFit[analModels.size()-1]){
            return analModels.size()-1; //use k-1 model
	 } else {
            return -1;  //no model selected
	 }
      }
   }
}

double calcSlopeFactor(double bmr, double bmdl){
   
   if (bmdl > 0.0){
      return bmr/bmdl;
   } else {
      return BMDS_MISSING;
   }
}

void BMDS_ENTRY_API __stdcall runMultitumorModel(struct python_multitumor_analysis *pyAnal, struct python_multitumor_result *pyRes){


   //calculate slopeFactor for all individual models
   for (int i=0; i<pyAnal->ndatasets; i++){
      for (int j=0; j<pyAnal->nmodels[i]; j++){
         double bmr = pyAnal->models[i][j].BMR;
         pyRes->models[i][j].bmdsRes.setSlopeFactor(bmr);
      }
   }


   //calculate combined LL & constant
   pyRes->combined_LL = 0.0; 
   pyRes->combined_LL_const = 0.0;
   for (int i=0; i<pyAnal->ndatasets; i++){
     int selIndex = pyRes->selectedModelIndex[i];
     pyRes->combined_LL += pyRes->models[i][selIndex].max;
     pyRes->combined_LL_const += LogLik_Constant(pyAnal->models[i][selIndex].Y, pyAnal->models[i][selIndex].n_group);
   }
   pyRes->combined_LL = pyRes->combined_LL*-1;
   pyRes->combined_LL_const = pyRes->combined_LL_const*-1;


   Multistage_ComboBMD(pyAnal, pyRes);

   pyRes->setSlopeFactor(pyAnal->BMR);

}





void BMDS_ENTRY_API __stdcall pythonBMDSNested(struct python_nested_analysis *pyAnal, struct python_nested_result *pyRes){
   
   //stubbed results
   pyRes->bmdsRes.validResult = true;
   pyRes->bmdsRes.BMD = 12.95166613;
   pyRes->bmdsRes.BMDL = 9.643478831;
   pyRes->bmdsRes.BMDU = (double)BMDS_MISSING;

   pyRes->bmdsRes.AIC = 546.9572409;
   pyRes->bmdsRes.chisq = 19.6087053;
   pyRes->combPVal = 0.994;
   pyRes->df = 35;

   pyRes->nparms = 9;
   pyRes->parms.push_back(0.084733516);
   pyRes->parms.push_back(-4.109651594);
   pyRes->parms.push_back(0.004761366);
   pyRes->parms.push_back(-0.055489253);
   pyRes->parms.push_back(1.0);
   pyRes->parms.push_back(0);
   pyRes->parms.push_back(0);
   pyRes->parms.push_back(0);
   pyRes->parms.push_back(0);

   pyAnal->seed = 1687267999;
   //pyAnal->iterations = 1000;
   pyRes->LL = -269.478205;
   pyRes->obsChiSq = 19.6087053;

   pyRes->boot.numRuns = 3;   
   pyRes->boot.pVal.push_back(0.994);
   pyRes->boot.pVal.push_back(0.997);
   pyRes->boot.pVal.push_back(0.991);
   pyRes->boot.pVal.push_back(0.994);
   pyRes->boot.perc50.push_back(38.79787451);
   pyRes->boot.perc50.push_back(38.19554318);
   pyRes->boot.perc50.push_back(37.75643018);
   pyRes->boot.perc50.push_back(38.26776644);
   pyRes->boot.perc90.push_back(51.1848655);
   pyRes->boot.perc90.push_back(50.4569082);
   pyRes->boot.perc90.push_back(50.3883043);
   pyRes->boot.perc90.push_back(50.6118335);
   pyRes->boot.perc95.push_back(54.69182);
   pyRes->boot.perc95.push_back(53.99859);
   pyRes->boot.perc95.push_back(54.18472);
   pyRes->boot.perc95.push_back(54.55846);
   pyRes->boot.perc99.push_back(60.479415);
   pyRes->boot.perc99.push_back(63.639965);
   pyRes->boot.perc99.push_back(61.778094);
   pyRes->boot.perc99.push_back(62.421371);

   pyRes->SRs.push_back(-0.31484);
   pyRes->SRs.push_back(0.314837);
   pyRes->SRs.push_back(-0.31484);
   pyRes->SRs.push_back(0.314837);
   pyRes->SRs.push_back(-0.31484);
   pyRes->SRs.push_back(0.314837);

   pyRes->litter.dose.push_back(0);
   pyRes->litter.dose.push_back(0);
   pyRes->litter.dose.push_back(0);
   pyRes->litter.dose.push_back(0);
   pyRes->litter.dose.push_back(0);
   pyRes->litter.dose.push_back(0);
   pyRes->litter.dose.push_back(0);
   pyRes->litter.dose.push_back(0);
   pyRes->litter.dose.push_back(0);
   pyRes->litter.dose.push_back(0);
   pyRes->litter.dose.push_back(25);
   pyRes->litter.dose.push_back(25);
   pyRes->litter.dose.push_back(25);
   pyRes->litter.dose.push_back(25);
   pyRes->litter.dose.push_back(25);
   pyRes->litter.dose.push_back(25);
   pyRes->litter.dose.push_back(25);
   pyRes->litter.dose.push_back(25);
   pyRes->litter.dose.push_back(25);
   pyRes->litter.dose.push_back(25);
   pyRes->litter.dose.push_back(50);
   pyRes->litter.dose.push_back(50);
   pyRes->litter.dose.push_back(50);
   pyRes->litter.dose.push_back(50);
   pyRes->litter.dose.push_back(50);
   pyRes->litter.dose.push_back(50);
   pyRes->litter.dose.push_back(50);
   pyRes->litter.dose.push_back(50);
   pyRes->litter.dose.push_back(50);
   pyRes->litter.dose.push_back(50);
   pyRes->litter.dose.push_back(100);
   pyRes->litter.dose.push_back(100);
   pyRes->litter.dose.push_back(100);
   pyRes->litter.dose.push_back(100);
   pyRes->litter.dose.push_back(100);
   pyRes->litter.dose.push_back(100);
   pyRes->litter.dose.push_back(100);
   pyRes->litter.dose.push_back(100);
   pyRes->litter.dose.push_back(100);

   pyRes->litter.LSC.push_back(9);
   pyRes->litter.LSC.push_back(9);
   pyRes->litter.LSC.push_back(10);
   pyRes->litter.LSC.push_back(10);
   pyRes->litter.LSC.push_back(11);
   pyRes->litter.LSC.push_back(13);
   pyRes->litter.LSC.push_back(14);
   pyRes->litter.LSC.push_back(14);
   pyRes->litter.LSC.push_back(15);
   pyRes->litter.LSC.push_back(16);
   pyRes->litter.LSC.push_back(9);
   pyRes->litter.LSC.push_back(9);
   pyRes->litter.LSC.push_back(10);
   pyRes->litter.LSC.push_back(10);
   pyRes->litter.LSC.push_back(11);
   pyRes->litter.LSC.push_back(12);
   pyRes->litter.LSC.push_back(13);
   pyRes->litter.LSC.push_back(14);
   pyRes->litter.LSC.push_back(14);
   pyRes->litter.LSC.push_back(14);
   pyRes->litter.LSC.push_back(7);
   pyRes->litter.LSC.push_back(10);
   pyRes->litter.LSC.push_back(10);
   pyRes->litter.LSC.push_back(11);
   pyRes->litter.LSC.push_back(11);
   pyRes->litter.LSC.push_back(11);
   pyRes->litter.LSC.push_back(11);
   pyRes->litter.LSC.push_back(14);
   pyRes->litter.LSC.push_back(14);
   pyRes->litter.LSC.push_back(15);
   pyRes->litter.LSC.push_back(8);
   pyRes->litter.LSC.push_back(10);
   pyRes->litter.LSC.push_back(11);
   pyRes->litter.LSC.push_back(11);
   pyRes->litter.LSC.push_back(12);
   pyRes->litter.LSC.push_back(12);
   pyRes->litter.LSC.push_back(13);
   pyRes->litter.LSC.push_back(14);
   pyRes->litter.LSC.push_back(14);

   pyRes->litter.estProb.push_back(0.127585814);
   pyRes->litter.estProb.push_back(0.127585814);
   pyRes->litter.estProb.push_back(0.13234718);
   pyRes->litter.estProb.push_back(0.13234718);
   pyRes->litter.estProb.push_back(0.137108547);
   pyRes->litter.estProb.push_back(0.14663128);
   pyRes->litter.estProb.push_back(0.151392646);
   pyRes->litter.estProb.push_back(0.151392646);
   pyRes->litter.estProb.push_back(0.156154013);
   pyRes->litter.estProb.push_back(0.160915379);
   pyRes->litter.estProb.push_back(0.301527034);
   pyRes->litter.estProb.push_back(0.301527034);
   pyRes->litter.estProb.push_back(0.297781775);
   pyRes->litter.estProb.push_back(0.297781775);
   pyRes->litter.estProb.push_back(0.294373053);
   pyRes->litter.estProb.push_back(0.291294677);
   pyRes->litter.estProb.push_back(0.288539894);
   pyRes->litter.estProb.push_back(0.286101463);
   pyRes->litter.estProb.push_back(0.286101463);
   pyRes->litter.estProb.push_back(0.286101463);
   pyRes->litter.estProb.push_back(0.433391608);
   pyRes->litter.estProb.push_back(0.410232266);
   pyRes->litter.estProb.push_back(0.410232266);
   pyRes->litter.estProb.push_back(0.40315061);
   pyRes->litter.estProb.push_back(0.40315061);
   pyRes->litter.estProb.push_back(0.40315061);
   pyRes->litter.estProb.push_back(0.40315061);
   pyRes->litter.estProb.push_back(0.38390157);
   pyRes->litter.estProb.push_back(0.38390157);
   pyRes->litter.estProb.push_back(0.378161814);
   pyRes->litter.estProb.push_back(0.572726279);
   pyRes->litter.estProb.push_back(0.553298381);
   pyRes->litter.estProb.push_back(0.543802864);
   pyRes->litter.estProb.push_back(0.543802864);
   pyRes->litter.estProb.push_back(0.534476942);
   pyRes->litter.estProb.push_back(0.534476942);
   pyRes->litter.estProb.push_back(0.52533783);
   pyRes->litter.estProb.push_back(0.516402011);
   pyRes->litter.estProb.push_back(0.516402011);

   pyRes->litter.litterSize.push_back(9);
   pyRes->litter.litterSize.push_back(9);
   pyRes->litter.litterSize.push_back(10);
   pyRes->litter.litterSize.push_back(10);
   pyRes->litter.litterSize.push_back(11);
   pyRes->litter.litterSize.push_back(13);
   pyRes->litter.litterSize.push_back(14);
   pyRes->litter.litterSize.push_back(14);
   pyRes->litter.litterSize.push_back(15);
   pyRes->litter.litterSize.push_back(16);
   pyRes->litter.litterSize.push_back(9);
   pyRes->litter.litterSize.push_back(9);
   pyRes->litter.litterSize.push_back(10);
   pyRes->litter.litterSize.push_back(10);
   pyRes->litter.litterSize.push_back(11);
   pyRes->litter.litterSize.push_back(12);
   pyRes->litter.litterSize.push_back(13);
   pyRes->litter.litterSize.push_back(14);
   pyRes->litter.litterSize.push_back(14);
   pyRes->litter.litterSize.push_back(14);
   pyRes->litter.litterSize.push_back(7);
   pyRes->litter.litterSize.push_back(10);
   pyRes->litter.litterSize.push_back(10);
   pyRes->litter.litterSize.push_back(11);
   pyRes->litter.litterSize.push_back(11);
   pyRes->litter.litterSize.push_back(11);
   pyRes->litter.litterSize.push_back(11);
   pyRes->litter.litterSize.push_back(14);
   pyRes->litter.litterSize.push_back(14);
   pyRes->litter.litterSize.push_back(15);
   pyRes->litter.litterSize.push_back(8);
   pyRes->litter.litterSize.push_back(10);
   pyRes->litter.litterSize.push_back(11);
   pyRes->litter.litterSize.push_back(11);
   pyRes->litter.litterSize.push_back(12);
   pyRes->litter.litterSize.push_back(12);
   pyRes->litter.litterSize.push_back(13);
   pyRes->litter.litterSize.push_back(14);
   pyRes->litter.litterSize.push_back(14);   

   pyRes->litter.expected.push_back(1.148272326);
   pyRes->litter.expected.push_back(1.148272326);
   pyRes->litter.expected.push_back(1.323471805);
   pyRes->litter.expected.push_back(1.323471805);
   pyRes->litter.expected.push_back(1.508194016);
   pyRes->litter.expected.push_back(1.906206638);
   pyRes->litter.expected.push_back(2.119497049);
   pyRes->litter.expected.push_back(2.119497049);
   pyRes->litter.expected.push_back(2.342310192);
   pyRes->litter.expected.push_back(2.574646069);
   pyRes->litter.expected.push_back(2.713743308);
   pyRes->litter.expected.push_back(2.713743308);
   pyRes->litter.expected.push_back(2.977817749);
   pyRes->litter.expected.push_back(2.977817749);
   pyRes->litter.expected.push_back(3.238103583);
   pyRes->litter.expected.push_back(3.495536119);
   pyRes->litter.expected.push_back(3.751018618);
   pyRes->litter.expected.push_back(4.005420479);
   pyRes->litter.expected.push_back(4.005420479);
   pyRes->litter.expected.push_back(4.005420479);
   pyRes->litter.expected.push_back(3.033741255);
   pyRes->litter.expected.push_back(4.102322662);
   pyRes->litter.expected.push_back(4.102322662);
   pyRes->litter.expected.push_back(4.434656714);
   pyRes->litter.expected.push_back(4.434656714);
   pyRes->litter.expected.push_back(4.434656714);
   pyRes->litter.expected.push_back(4.434656714);
   pyRes->litter.expected.push_back(5.374621982);
   pyRes->litter.expected.push_back(5.374621982);
   pyRes->litter.expected.push_back(5.672427209);
   pyRes->litter.expected.push_back(4.581810233);
   pyRes->litter.expected.push_back(5.532983808);
   pyRes->litter.expected.push_back(5.981831502);
   pyRes->litter.expected.push_back(5.981831502);
   pyRes->litter.expected.push_back(6.413723301);
   pyRes->litter.expected.push_back(6.413723301);
   pyRes->litter.expected.push_back(6.829391788);
   pyRes->litter.expected.push_back(7.229628148);
   pyRes->litter.expected.push_back(7.229628148);

   pyRes->litter.observed.push_back(0);
   pyRes->litter.observed.push_back(1);
   pyRes->litter.observed.push_back(1);
   pyRes->litter.observed.push_back(2);
   pyRes->litter.observed.push_back(2);
   pyRes->litter.observed.push_back(3);
   pyRes->litter.observed.push_back(3);
   pyRes->litter.observed.push_back(2);
   pyRes->litter.observed.push_back(2);
   pyRes->litter.observed.push_back(1);
   pyRes->litter.observed.push_back(2);
   pyRes->litter.observed.push_back(5);
   pyRes->litter.observed.push_back(2);
   pyRes->litter.observed.push_back(1);
   pyRes->litter.observed.push_back(4);
   pyRes->litter.observed.push_back(3);
   pyRes->litter.observed.push_back(6);
   pyRes->litter.observed.push_back(6);
   pyRes->litter.observed.push_back(4);
   pyRes->litter.observed.push_back(3);
   pyRes->litter.observed.push_back(2);
   pyRes->litter.observed.push_back(5);
   pyRes->litter.observed.push_back(5);
   pyRes->litter.observed.push_back(4);
   pyRes->litter.observed.push_back(4);
   pyRes->litter.observed.push_back(5);
   pyRes->litter.observed.push_back(4);
   pyRes->litter.observed.push_back(4);
   pyRes->litter.observed.push_back(5);
   pyRes->litter.observed.push_back(6);
   pyRes->litter.observed.push_back(5);
   pyRes->litter.observed.push_back(4);
   pyRes->litter.observed.push_back(6);
   pyRes->litter.observed.push_back(6);
   pyRes->litter.observed.push_back(8);
   pyRes->litter.observed.push_back(8);
   pyRes->litter.observed.push_back(7);
   pyRes->litter.observed.push_back(6);
   pyRes->litter.observed.push_back(6);

   pyRes->litter.SR.push_back(-1.147257987);
   pyRes->litter.SR.push_back(-0.148141348);
   pyRes->litter.SR.push_back(-0.301860366);
   pyRes->litter.SR.push_back(0.631328744);
   pyRes->litter.SR.push_back(0.431109028);
   pyRes->litter.SR.push_back(0.857594396);
   pyRes->litter.SR.push_back(0.65653973);
   pyRes->litter.SR.push_back(-0.089101984);
   pyRes->litter.SR.push_back(-0.243481535);
   pyRes->litter.SR.push_back(-1.071325161);
   pyRes->litter.SR.push_back(-0.518421334);
   pyRes->litter.SR.push_back(1.660602955);
   pyRes->litter.SR.push_back(-0.676196332);
   pyRes->litter.SR.push_back(-1.367732493);
   pyRes->litter.SR.push_back(0.504037656);
   pyRes->litter.SR.push_back(-0.314836859);
   pyRes->litter.SR.push_back(1.376689417);
   pyRes->litter.SR.push_back(1.179530167);
   pyRes->litter.SR.push_back(-0.003205497);
   pyRes->litter.SR.push_back(-0.594573329);
   pyRes->litter.SR.push_back(-0.788462565);
   pyRes->litter.SR.push_back(0.577118305);
   pyRes->litter.SR.push_back(0.577118305);
   pyRes->litter.SR.push_back(-0.267167737);
   pyRes->litter.SR.push_back(-0.267167737);
   pyRes->litter.SR.push_back(0.347496039);
   pyRes->litter.SR.push_back(-0.267167737);
   pyRes->litter.SR.push_back(-0.755412682);
   pyRes->litter.SR.push_back(-0.205870559);
   pyRes->litter.SR.push_back(0.174415333);
   pyRes->litter.SR.push_back(0.298883377);
   pyRes->litter.SR.push_back(-0.975099884);
   pyRes->litter.SR.push_back(0.010998302);
   pyRes->litter.SR.push_back(0.010998302);
   pyRes->litter.SR.push_back(0.918022312);
   pyRes->litter.SR.push_back(0.918022312);
   pyRes->litter.SR.push_back(0.094758157);
   pyRes->litter.SR.push_back(-0.65761782);
   pyRes->litter.SR.push_back(-0.65761782);

   pyRes->reduced.dose.push_back(0);
   pyRes->reduced.dose.push_back(25);
   pyRes->reduced.dose.push_back(50);
   pyRes->reduced.dose.push_back(100);
   pyRes->reduced.propAffect.push_back(0.14049587);
   pyRes->reduced.propAffect.push_back(0.31034482);
   pyRes->reduced.propAffect.push_back(0.38596491);
   pyRes->reduced.propAffect.push_back(0.53333333);
   pyRes->reduced.lowerConf.push_back(0.10100315);
   pyRes->reduced.lowerConf.push_back(0.23266263);
   pyRes->reduced.lowerConf.push_back(0.34156249);
   pyRes->reduced.lowerConf.push_back(0.46300766);
   pyRes->reduced.upperConf.push_back(0.19144442);
   pyRes->reduced.upperConf.push_back(0.39977302);
   pyRes->reduced.upperConf.push_back(0.43230026);
   pyRes->reduced.upperConf.push_back(0.60240216);

}


double round_to(double value, double precision ){
   double res = std::round(value / precision) * precision;
   return res;
}


