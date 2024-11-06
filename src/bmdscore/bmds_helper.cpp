#ifdef WIN32
#  include "pch.h"
#else
#  include "stdafx.h"
#endif
#include <math.h>
#include <stdio.h>

#include <ctime>
#include <iomanip>
// #include <cmath>
#include <nlopt.hpp>

#include "analysis_of_deviance.h"
#include "bmds_helper.h"

// calendar versioning; see https://peps.python.org/pep-0440/#pre-releases
std::string BMDS_VERSION = "24.1a4";

double python_dichotomous_model_result::getSRAtDose(double targetDose, std::vector<double> doses) {
  std::vector<double> diff;
  double absDiff = DBL_MAX;
  double srVal = BMDS_MISSING;
  if (!bmdsRes.validResult || doses.size() != gof.residual.size()) {
    return BMDS_MISSING;
  }
  for (int i = 0; i < doses.size(); i++) {
    diff.push_back(abs(targetDose - doses[i]));
  }
  int minIndex =
      std::distance(std::begin(diff), std::min_element(std::begin(diff), std::end(diff)));

  return gof.residual[minIndex];
}

int checkForBoundedParms(
    int nparms, double *parms, double *lowerBound, double *upperBound, struct BMDS_results *BMDSres
) {
  // First find number of bounded parms
  int bounded = 0;
  BMDSres->bounded.clear();
  for (int i = 0; i < nparms; i++) {
    BMDSres->bounded.push_back(false);
    // 5*i+4 is location of min in prior array
    // 5*i+5 is location of max in prior array
    if (fabs(parms[i] - lowerBound[i]) < BMDS_EPS || fabs(parms[i] - upperBound[i]) < BMDS_EPS) {
      bounded++;
      BMDSres->bounded[i] = true;
    }
  }
  return bounded;
}

void calcDichoAIC(
    struct dichotomous_analysis *anal, struct dichotomous_model_result *res,
    struct BMDS_results *BMDSres, double estParmCount
) {
  bool freqModel = anal->prior[0] == 0;

  //  // First find number of bounded parms
  // int bounded = checkForBoundedParms(anal->parms, res->parms, anal->prior, BMDSres);

  double *lowerBound = (double *)malloc(anal->parms * sizeof(double));
  double *upperBound = (double *)malloc(anal->parms * sizeof(double));

  // copy prior values
  for (int i = 0; i < anal->parms; i++) {
    lowerBound[i] = anal->prior[3 * anal->parms + i];
    upperBound[i] = anal->prior[4 * anal->parms + i];
  }

  // scale priors as needed
  rescale_dichoParms(anal->model, lowerBound);
  rescale_dichoParms(anal->model, upperBound);

  int bounded = checkForBoundedParms(anal->parms, res->parms, lowerBound, upperBound, BMDSres);

  estParmCount = anal->parms - bounded;

  // if freq then model_df should be rounded to nearest whole number
  if (freqModel) estParmCount = round(estParmCount);
  BMDSres->AIC = 2 * (res->max + estParmCount);

  free(lowerBound);
  free(upperBound);
}

void calcContAIC(
    struct continuous_analysis *anal, struct continuous_model_result *res,
    struct BMDS_results *BMDSres
) {
  bool freqModel = anal->prior[0] == 0;

  // First find number of bounded parms
  double *lowerBound = (double *)malloc(anal->parms * sizeof(double));
  double *upperBound = (double *)malloc(anal->parms * sizeof(double));

  // copy prior values
  for (int i = 0; i < anal->parms; i++) {
    lowerBound[i] = anal->prior[3 * anal->parms + i];
    upperBound[i] = anal->prior[4 * anal->parms + i];
  }
  // scale priors as needed
  rescale_contParms(anal, lowerBound);
  rescale_contParms(anal, upperBound);

  if (anal->model == cont_model::exp_3) {
    for (int i = 2; i < anal->parms - 1; i++) {
      lowerBound[i] = lowerBound[i + 1];
      upperBound[i] = upperBound[i + 1];
    }
  }
  int bounded = checkForBoundedParms(res->nparms, res->parms, lowerBound, upperBound, BMDSres);
  // int bounded = checkForBoundedParms(anal->parms, res->parms, lowerBound, upperBound, BMDSres);

  double estParmCount = res->model_df - bounded;
  // if freq then model_df should be rounded to nearest whole number
  if (freqModel) estParmCount = round(estParmCount);

  BMDSres->AIC = 2 * (res->max + estParmCount);

  //               //    BMDSres->AIC = -9998.0;
  free(lowerBound);
  free(upperBound);
}

double findQuantileVals(double *quant, double *val, int arrSize, double target) {
  double retVal = BMDS_MISSING;

  for (int i = 0; i < arrSize; i++) {
    if (fabs(quant[i] - target) < BMDS_EPS && std::isfinite(val[i])) {
      // exact match
      retVal = val[i];
      break;
    } else if (quant[i] > target && i > 0) {
      // linear interpolation
      retVal = val[i - 1] +
               ((val[i] - val[i - 1]) / (quant[i] - quant[i - 1])) * (target - quant[i - 1]);
      break;
    }
  }
  return retVal;
}

void collect_dicho_bmd_values(
    struct dichotomous_analysis *anal, struct dichotomous_model_result *res,
    struct BMDS_results *BMDSres, double estParmCount
) {
  int distSize = res->dist_numE * 2;

  double *quant = (double *)malloc(distSize / 2 * sizeof(double));
  double *val = (double *)malloc(distSize / 2 * sizeof(double));

  for (int i = 0; i < distSize / 2; i++) {
    val[i] = res->bmd_dist[i];
  }
  for (int i = distSize / 2; i < distSize; i++) {
    quant[i - distSize / 2] = res->bmd_dist[i];
  }

  calcDichoAIC(anal, res, BMDSres, estParmCount);
  BMDSres->BMD = findQuantileVals(quant, val, distSize / 2, 0.50);
  BMDSres->BMDL = findQuantileVals(quant, val, distSize / 2, anal->alpha);
  BMDSres->BMDU = findQuantileVals(quant, val, distSize / 2, 1.0 - anal->alpha);

  free(quant);
  free(val);
}

void collect_dichoMA_bmd_values(
    struct dichotomousMA_analysis *anal, struct dichotomousMA_result *res,
    struct BMDSMA_results *BMDSres, double alpha
) {
  int distSize = res->dist_numE * 2;

  double *quantMA = (double *)malloc(distSize / 2 * sizeof(double));
  double *valMA = (double *)malloc(distSize / 2 * sizeof(double));

  for (int i = 0; i < distSize / 2; i++) {
    valMA[i] = res->bmd_dist[i];
  }
  for (int i = distSize / 2; i < distSize; i++) {
    quantMA[i - distSize / 2] = res->bmd_dist[i];
  }

  //  calculate MA quantiles
  BMDSres->BMD_MA = findQuantileVals(quantMA, valMA, distSize / 2, 0.50);
  BMDSres->BMDL_MA = findQuantileVals(quantMA, valMA, distSize / 2, alpha);
  BMDSres->BMDU_MA = findQuantileVals(quantMA, valMA, distSize / 2, 1.0 - alpha);

  // calculate individual model quantiles
  for (int j = 0; j < anal->nmodels; j++) {
    for (int i = 0; i < distSize / 2; i++) {
      valMA[i] = res->models[j]->bmd_dist[i];
    }
    for (int i = distSize / 2; i < distSize; i++) {
      quantMA[i - distSize / 2] = res->models[j]->bmd_dist[i];
    }
    BMDSres->BMD[j] = findQuantileVals(quantMA, valMA, distSize / 2, 0.50);
    BMDSres->BMDL[j] = findQuantileVals(quantMA, valMA, distSize / 2, alpha);
    BMDSres->BMDU[j] = findQuantileVals(quantMA, valMA, distSize / 2, 1.0 - alpha);
  }
  free(quantMA);
  free(valMA);
}

void collect_cont_bmd_values(
    struct continuous_analysis *anal, struct continuous_model_result *res,
    struct BMDS_results *BMDSres
) {
  int distSize = res->dist_numE * 2;

  double *contVal = (double *)malloc(distSize / 2 * sizeof(double));
  double *contQuant = (double *)malloc(distSize / 2 * sizeof(double));

  for (int i = 0; i < distSize / 2; i++) {
    contVal[i] = res->bmd_dist[i];
  }

  for (int i = distSize / 2; i < distSize; i++) {
    contQuant[i - distSize / 2] = res->bmd_dist[i];
  }

  calcContAIC(anal, res, BMDSres);
  BMDSres->BMD = findQuantileVals(contQuant, contVal, distSize / 2, 0.50);
  BMDSres->BMDL = findQuantileVals(contQuant, contVal, distSize / 2, anal->alpha);
  BMDSres->BMDU = findQuantileVals(contQuant, contVal, distSize / 2, 1.0 - anal->alpha);

  free(contVal);
  free(contQuant);
}

/*****************************************************************
 * DLgama - double log gamma
 *  input:
 *   X - double value
 *  Returns log(z-1)!
 *************************/
double DLgamma(double x) {
  double r, az, z, y, dl;
  if (x == 0) return 1;
  z = x;
  r = 12.0;
  az = 1.0;
  for (; z - r < 0; z++) az = az * z;
  y = z * z;
  dl = 0.91893853320467274 - z + (z - 0.5) * log(z) +
       z * (0.08333333333333333 * y + 0.0648322851140734) /
           (y * (y + 0.811320754825416) + 0.01752021563249146);
  dl = dl - log(az);
  return dl;
}

// calculate dichotomous log likelihood constant
double LogLik_Constant(std::vector<double> Y, std::vector<double> n_group) {
  double X, N, NX, con, xf, nf, nxf, val;

  con = 0;
  for (int i = 0; i < Y.size(); i++) {
    X = Y[i] + 1;
    N = n_group[i] + 1;
    NX = n_group[i] - Y[i] + 1;
    xf = exp(DLgamma(X));
    nf = exp(DLgamma(N));
    nxf = exp(DLgamma(NX));
    val = log(nf / (xf * nxf));
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

double zeroin(
    double ax, double bx, double tol, double (*f)(int, double[], double, double), int nparm,
    double Parms[], double ck
)
/* ax        Left border | of the range */
/* bx        Right border | the root is sought*/
/* f	  Function under investigation */
/* nparm     number of parameters in Parms */
/* Parms     vector of parameters to pass to f */
/* tol       Acceptable tolerance for the root */
/* gtol      tolerance to pass to f */
// TODO: clean up comments after debugging is complete
{
  double a, b, c; /* Abscissae, descr. see above	*/
  double fa;      /* f(a)				*/
  double fb;      /* f(b)				*/
  double fc;      /* f(c)				*/

  std::cout << std::fixed << std::showpoint;
  std::cout << std::setprecision(15);

  a = ax;
  b = bx;
  fa = (*f)(nparm - 1, Parms, a, ck);
  fb = (*f)(nparm - 1, Parms, b, ck);
  c = a;
  fc = fa;

  for (;;) /* Main iteration loop	*/
  {
    double prev_step = b - a; /* Distance from the last but one*/
                              /* to the last approximation	*/
    double tol_act;           /* Actual tolerance		*/
    double p;                 /* Interpolation step is calcu- */
    double q;                 /* lated in the form p/q; divi- */
                              /* sion operations is delayed   */
    /* until the last moment	*/
    double new_step;           /* Step at this iteration       */
    if (fabs(fc) < fabs(fb)) { /* Swap data for b to be the 	*/
      a = b;
      b = c;
      c = a; /* best approximation		*/
      fa = fb;
      fb = fc;
      fc = fa;
    }
    tol_act = 2 * DBL_EPSILON * fabs(b) + tol / 2;
    new_step = (c - b) / 2;

    if (fabs(new_step) <= tol_act || fb == (double)0) {
      return b; /* Acceptable approx. is found	*/
    }

    /* Decide if the interpolation can be tried	*/
    if (fabs(prev_step) >= tol_act /* If prev_step was large enough*/
        && fabs(fa) > fabs(fb))    /* and was in true direction,	*/
    {                              /* Interpolatiom may be tried	*/
      double t1, cb, t2;
      cb = c - b;
      if (a == c)     /* If we have only two distinct	*/
      {               /* points linear interpolation 	*/
        t1 = fb / fa; /* can only be applied		*/
        p = cb * t1;
        q = 1.0 - t1;
      } else /* Quadric inverse interpolation*/
      {
        q = fa / fc;
        t1 = fb / fc;
        t2 = fb / fa;
        p = t2 * (cb * q * (q - t1) - (b - a) * (t1 - 1.0));
        q = (q - 1.0) * (t1 - 1.0) * (t2 - 1.0);
      }
      if (p > (double)0) /* p was calculated with the op-*/
        q = -q;          /* posite sign; make p positive	*/
      else               /* and assign possible minus to	*/
        p = -p;          /* q				*/

      if (p < (0.75 * cb * q - fabs(tol_act * q) / 2) /* If b+p/q falls in [b,c]*/
          && p < fabs(prev_step * q / 2))             /* and isn't too large	*/
        new_step = p / q;                             /* it is accepted	*/
                                                      /* If p/q is too large then the	*/
                                                      /* bissection procedure can 	*/
                                                      /* reduce [b,c] range to more	*/
                                                      /* extent			*/
    }

    if (fabs(new_step) < tol_act) /* Adjust the step to be not less*/
    {
      if (new_step > (double)0) /* than tolerance		*/
        new_step = tol_act;
      else
        new_step = -tol_act;
    }

    a = b;
    fa = fb; /* Save the previous approx.	*/
    b += new_step;
    fb = (*f)(nparm - 1, Parms, b, ck);             /* Do step to a new approxim.	*/
    if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) { /* Adjust c for it to have a sign*/
      c = a;
      fc = fa; /* opposite to that of b	*/
    }
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
double BMD_func(int n, double p[], double x, double ck) {
  double poly, fx, D;
  int j;
  D = exp(x);
  poly = p[n];
  for (j = n - 1; j > 0; j--) {
    poly = poly * D + p[j];
  }
  poly = poly * D;
  fx = poly - ck; /* ck = log(1-A) */
  return fx;
}

double getclmt(
    python_multitumor_analysis *pyAnal, python_multitumor_result *pyRes, double Dose, double target,
    double maxDose, std::vector<double> xParms, bool isBMDL
) {
  int nT = pyRes->selectedModelIndex.size();
  double bmr = pyAnal->BMR;
  std::vector<int> degree;
  for (int j = 0; j < nT; j++) {
    int modDeg = pyAnal->models[j][pyRes->selectedModelIndex[j]].degree;
    degree.push_back(modDeg);
  }
  int N = xParms.size() + 1;
  double bmd = Dose;
  std::vector<double> x(N);
  int iOffset = 0;
  x[0] = log(bmd);

  // get max degree for loop
  int maxDegree = *max_element(degree.begin(), degree.end());

  std::vector<std::vector<double>> tmp2;
  // restructure to group like terms together
  int count = 0;
  for (int i = 0; i < nT; i++) {
    std::vector<double> theParms;
    for (int j = 0; j <= degree[i]; j++) {
      theParms.push_back(xParms[count]);
      count++;
    }
    tmp2.push_back(theParms);
  }

  count = 1;
  for (int j = 0; j <= maxDegree; j++) {
    for (int i = 0; i < nT; i++) {
      if (j < tmp2[i].size()) {
        x[count] = tmp2[i][j];
        count++;
      }
    }
  }

  // need to round to roughly single-precision for convergence?????
  for (int j = 1; j < x.size(); j++) {
    x[j] = round_to(x[j], 0.00001);
  }

  std::vector<double> lb(x.size());
  std::vector<double> ub(x.size());

  // calculate the values that will correspond to '0' and 'infinity' in BMD calcs
  double lminbmd = log(DBL_MIN) - log(maxDose);
  double lmaxbmd = log(DBL_MAX) - log(maxDose);

  lb[0] = lminbmd;  // BMD lower limit
  for (int i = 1; i < x.size(); i++) {
    lb[i] = 0.0;  // beta min value
  }

  ub[0] = log(maxDose);
  for (int i = 1; i < x.size(); i++) {
    ub[i] = 1e4;  // beta max value
  }

  //   //constraint data
  struct msComboEq eq1;
  eq1.bmr = bmr;
  eq1.nT = nT;
  eq1.degree = degree;

  struct msComboInEq ineq1;
  ineq1.nT = nT;
  ineq1.target = target;
  // TODO need to add handling of failed datasets

  for (int i = 0; i < pyRes->selectedModelIndex.size(); i++) {
    int selIndex = pyRes->selectedModelIndex[i];
    std::vector<double> scaledDose = pyAnal->models[i][selIndex].doses;
    for (int j = 0; j < scaledDose.size(); j++) {
      scaledDose[j] /= maxDose;
    }
    ineq1.doses.push_back(scaledDose);
    ineq1.Y.push_back(pyAnal->models[i][selIndex].Y);
    ineq1.n_group.push_back(pyAnal->models[i][selIndex].n_group);
  }
  ineq1.nObs = pyAnal->n;
  ineq1.degree = degree;

  nlopt::opt opt(nlopt::LD_SLSQP, x.size());
  if (isBMDL) {
    opt.set_min_objective(objfunc_bmdl, NULL);
  } else {
    opt.set_min_objective(objfunc_bmdu, NULL);
  }
  opt.add_equality_constraint(myEqualityConstraint, &eq1, 1e-8);
  opt.add_inequality_constraint(myInequalityConstraint1, &ineq1, 1e-8);
  opt.set_xtol_rel(1e-8);
  opt.set_maxeval(100000);
  opt.set_lower_bounds(lb);
  opt.set_upper_bounds(ub);

  double minf, val;
  nlopt::result result = nlopt::FAILURE;
  try {
    result = opt.optimize(x, minf);
    //     std::cout << "found minimum at f(" << x[0] << ") = " << std::setprecision(10) << minf <<
    //     std::endl;
  } catch (std::exception &e) {
    std::cout << "nlopt failed: " << e.what() << std::endl;
  }

  val = x[0];

  return val;
}

/*****************************************************************
 * multitumorCLs -- returns the lower & upper confidence limits (BMDL & BMDU).
 *****************************************************************/
// TODO: combine with BMDU_combofunc
void multitumorCLs(
    struct python_multitumor_analysis *pyAnal, struct python_multitumor_result *pyRes, double Dose,
    double D, double LR, double gtol, int *is_zero
) { /* ck and LR are calculated in Multistage_ComboBMD() */

  int lnParmMax;
  double bmdl, bmdu, target, xmax;
  int i, j, nCall, k, nParmMax, ii = 0;
  double scale;

  std::vector<double> adxmax;
  int nParms;
  bmdl = bmdu = target = xmax = 0.0;

  adxmax.resize(pyRes->ndatasets, 0.0);

  lnParmMax = 0;
  nParms = 0;
  for (i = 0; i < pyRes->ndatasets; i++) {
    int selModelIndex = pyRes->selectedModelIndex[i];
    struct python_dichotomous_analysis mod = pyAnal->models[i][selModelIndex];
    struct python_dichotomous_model_result modRes = pyRes->models[i][selModelIndex];
    xmax = mod.doses[0];
    nParms += modRes.nparms;

    if (modRes.nparms > lnParmMax) {
      lnParmMax = modRes.nparms;
    }

    for (j = 0; j < mod.n; j++) {
      if (mod.doses[j] > xmax) {
        xmax = mod.doses[j];
      }
    }
    adxmax[i] = xmax;
  }

  /** rescale all doses to be: 0 <= Dose <= 1 **/
  xmax = adxmax[0];
  for (i = 0; i < pyRes->ndatasets; i++) {
    if (adxmax[i] > xmax) xmax = adxmax[i];
  }
  scale = xmax;

  nParmMax = (int)lnParmMax;
  std::vector<double> pdParms(nParms, 0.0);
  std::vector<double> pdParmsBak(nParms, 0.0);
  std::vector<double> pdParms2(nParms, 0.0);
  std::vector<double> pdVals(nParms, 0.0);
  std::vector<int> piSpec2(nParms, 0.0);

  k = -1;
  for (j = 0; j < nParmMax; j++) {
    for (i = pyAnal->ndatasets - 1; i >= 0; i--) {
      if (j < pyRes->models[i][0].nparms) {
        k++;
        pdParms[k] = pyRes->models[i][0].parms[j];
      }
    }
  }

  i = 0;
  j = 0;

  Dose = Dose / scale;

  target = (pyRes->combined_LL - LR); /* The value we want the likelihood */

  k = -1;

  for (j = 0; j < nParmMax; j++) {
    for (i = pyRes->ndatasets - 1; i >= 0; i--) {
      int iParms = pyRes->models[i][0].nparms;

      if (j < iParms) {
        k++;
        pdParmsBak[k] = pyRes->models[i][0].parms[j];
        pdParms[k] = pyRes->models[i][0].parms[j] * (pow(scale, (j)));
      }
    }
  }
  //  /* One more step for the background terms */
  for (i = 0; i < pyRes->ndatasets; i++) {
    pdParms[i] = -log(1.0 - pdParms[i]);
  }
  i = 0;

  bmdl = getclmt(pyAnal, pyRes, Dose, target, xmax, pdParms, true);

  bmdu = getclmt(pyAnal, pyRes, Dose, target, xmax, pdParms, false);

  pdParms2 = pdParms;
  int nT = 0;
  for (int i = 0; i < pyRes->ndatasets; i++) {
    nT++;
  }
  std::vector<std::vector<double>> ppdParms(nT, std::vector<double>(nParmMax, BMDS_MISSING));

  k = -1;
  for (int j = 0; j < nParmMax; j++) {
    for (int i = 0; i < pyRes->ndatasets; i++) {
      if (j < pyRes->models[i][0].nparms) {
        k++;
        ppdParms[i][j] = pdParms2[k];
        if (k < pyAnal->ndatasets) {
          ppdParms[i][j] = 1 - exp(-pdParms2[k]);
        }
      }
    }
  }

  int flag = 1;
  bmdl = exp(bmdl) * scale;
  bmdu = exp(bmdu) * scale;

  pyRes->BMDL = bmdl;
  pyRes->BMDU = bmdu;
  // xlk3 = ComboMaxLike2(flag,bmdl,&crisk, ppdParms, pyAnal, pyRes);
}

double ComboMaxLike2(
    int flag, double dose, double *crisk, std::vector<std::vector<double>> p,
    python_multitumor_analysis *pyAnal, python_multitumor_result *pyRes
) {
  int nObs, nT, nParms;
  double prob, like, dSumParm1, pr, bkg;

  prob = like = dSumParm1 = 0.0;

  for (int n = 0; n < pyRes->ndatasets; n++) {
    dSumParm1 += p[n][0];
    nObs = pyAnal->models[n][0].n;
    int selIndex = pyRes->selectedModelIndex[n];
    nParms = pyRes->models[n][selIndex].nparms;
    for (int i = 0; i < nObs; i++) {
      double D, Yn, Yp;
      D = pyAnal->models[n][selIndex].doses[i];
      Yp = pyAnal->models[n][selIndex].Y[i];
      Yn = pyAnal->models[n][selIndex].n_group[i] - Yp;
      prob = p[n][nParms - 1];
      if (n == 0 && flag == 1 && nParms == 2) {
        prob = p[0][1];
        for (int nt = 1; nt < pyRes->ndatasets; nt++) {
          prob -= p[nt][1];
        }
      }

      for (int j = nParms - 1; j >= 0; j--) {
        if (n == 0 && flag == 1 && j == 1) {
          pr = p[0][1];
          for (int nt = 1; nt < pyRes->ndatasets; nt++) {
            pr -= p[nt][1];
          }
          prob = D * prob + pr;
        } else {
          prob = D * prob + p[n][j];
        }
      }
      prob = (1 - exp(-1.0 * prob));
      if ((prob == 0) || (prob == 1)) {
        if (Yp <= 0 || Yn <= 0) {
          like += 0;
        } else {
          if (prob == 1) {
            like += Yn * (-710);
          } else {
            like += Yp * (-710);
          }
        }
      } else {
        like += Yp * log(prob) + Yn * log(1 - prob);
      }
    }
  }

  bkg = 1.0 - exp(-1.0 * (dSumParm1));

  for (int n = 0; n < pyRes->ndatasets; n++) {
    int selIndex = pyRes->selectedModelIndex[n];
    nParms = pyRes->models[n][selIndex].nparms;

    prob = p[n][nParms - 1];
    if (n == 0 && flag == 1 && nParms == 2) {
      prob = p[0][1];
      for (int nt = 1; nt < pyRes->ndatasets; nt++) {
        prob -= p[nt][1];
      }
    }
    for (int j = nParms - 1; j >= 0; j--) {
      if (n == 0 && flag == 1 && j == 1) {
        pr = p[0][1];
        for (int nt = 1; nt < pyRes->ndatasets; nt++) {
          pr -= p[nt][1];
        }
        prob = dose * prob + pr;
      } else {
        prob = dose * prob + p[n][j];
      }
    }
  }

  if (bkg == 1.0) {
    *crisk = 0.0;
  } else {
    *crisk = ((1.0 - exp(-1.0 * (prob))) - bkg) / (1.0 - bkg);
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
void Multistage_ComboBMD(
    struct python_multitumor_analysis *pyAnal, struct python_multitumor_result *pyRes
) {
  int cnparm, selIndex, nT, is_zero;
  double LR, xa, xb, D, fa, fb, Drange, cxmax, ck, poly;
  double gtol = 1e-12;
  double tol;  // for zeroin function
  std::vector<double> cp;

  nT = pyRes->ndatasets;
  // find largest nparm and largest dose
  cnparm = 0;
  cxmax = 0;
  for (int i = 0; i < nT; i++) {
    selIndex = pyRes->selectedModelIndex[i];
    if (pyRes->models[i][selIndex].nparms > cnparm) {
      cnparm = pyRes->models[i][selIndex].nparms;
    }
    double tmpMax = *max_element(
        std::begin(pyAnal->models[i][selIndex].doses), std::end(pyAnal->models[i][selIndex].doses)
    );
    if (cxmax < tmpMax) cxmax = tmpMax;
  }
  // add all model parameters to combined p
  cp.resize(cnparm, 0.0);
  for (int i = 0; i < nT; i++) {
    selIndex = pyRes->selectedModelIndex[i];
    cp[0] = cp[0] - log(1.0 - pyRes->models[i][selIndex].parms[0]);
    for (int j = 1; j < pyRes->models[i][selIndex].nparms; j++) {
      cp[j] = cp[j] + pyRes->models[i][selIndex].parms[j];
    }
  }
  // compute chi-squared value
  double cl = 1.0 - pyAnal->alpha;
  if (cl < 0.5) {
    LR = 0.5 * gsl_cdf_chisq_Pinv(1.0 - 2.0 * cl, 1);
  } else {
    LR = 0.5 * gsl_cdf_chisq_Pinv(2.0 * cl - 1.0, 1);
  }

  ck = -log(1 - pyAnal->BMR);  // Extra risk
  // solve the BMD
  xa = D = 0.0;
  fa = -ck;  // note:  ck>0.0
  fb = fa;
  Drange = cxmax;

  int k = 1;
  while (k < 300 && fb < 0) {
    fa = fb;
    xa = D;
    D = Drange * k / 100.0;
    poly = cp[cnparm - 1];
    for (int j = cnparm - 2; j > 0; j--) {
      poly = poly * D + cp[j];
    }
    poly = poly * D;
    fb = poly - ck;
    k++;
  }

  if (fb < 0)
    std::cout << "BMD Computation failed.  BMD is larger than three times maximum input doses."
              << std::endl;
  xb = D;
  tol = 1.0e-15;

  // compute BMD
  // BMD_func works on log scale, so convert xa and xb to logs
  if (xa == 0.0)
    xa = -690.7755;
  else
    xa = log(xa);
  xb = log(xb);

  xb = zeroin(xa, xb, tol, BMD_func, cnparm, &cp[0], ck);
  xa = exp(xa);
  xb = exp(xb);
  pyRes->BMD = xb;

  is_zero = 0;

  multitumorCLs(pyAnal, pyRes, xb, xa, LR, tol, &is_zero);
}

double objfunc_bmdl(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data) {
  // obj function of form F(BMD, beta) = BMD, where BMD=X[0] and beta=X[i] where i>0
  if (!grad.empty()) {
    // fill all gradients to zero (grad[1] to grad[n] are all zero
    //  because objective function only depends on X[0]
    std::fill(grad.begin(), grad.end(), 0);
    // set first grad to 1, since x[1] should be the BMD (Dose)
    grad[0] = 1.0;
  }
  return x[0];
}

double objfunc_bmdu(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data) {
  // obj function of form F(BMD, beta) = BMD, where BMD=X[0] and beta=X[i] where i>0
  if (!grad.empty()) {
    // fill all gradients to zero (grad[1] to grad[n] are all zero
    //  because objective function only depends on X[0]
    std::fill(grad.begin(), grad.end(), 0);
    // set first grad to 1, since x[1] should be the BMD (Dose)
    grad[0] = -1.0;
  }
  return -1 * x[0];
}

double myEqualityConstraint(const std::vector<double> &x, std::vector<double> &grad, void *data) {
  msComboEq *d = reinterpret_cast<msComboEq *>(data);
  double bmr = d->bmr;
  int nT = d->nT;
  std::vector<int> degree = d->degree;

  double D = exp(x[0]);
  int iIndex = x.size() - 1;
  double sum2 = 0.0;
  double sum = 0.0;

  if (!grad.empty()) {
    for (int l = nT - 1; l >= 0; l--) {
      sum = 0.0;
      for (int k = degree[l]; k > 0; k--) {
        sum = sum * D + k * x[iIndex] * D;
        iIndex -= 1;
      }
      iIndex -= 1;
      sum2 += sum;
    }
    grad[0] = sum2;

    iIndex = 1;
    for (int k = 0; k < nT; k++) {
      iIndex++;
      grad[iIndex] = D;
      for (int j = 2; j <= degree[k]; j++) {
        iIndex++;
        grad[iIndex] = pow(D, j);
      }
      iIndex++;
    }
  }

  // equality constraint calc
  sum = log(1.0 - bmr);
  iIndex = x.size() - 1;
  double sum3 = 0.0;
  for (int l = nT - 1; l >= 0; l--) {
    sum2 = 0.0;
    for (int k = degree[l]; k > 0; k--) {
      sum2 = sum2 * D + x[iIndex] * D;
      iIndex -= 1;
    }
    iIndex -= 1;
    sum3 += sum2;
  }

  return sum + sum3;
}

double myInequalityConstraint1(
    const std::vector<double> &x, std::vector<double> &grad, void *data
) {
  msComboInEq *d = reinterpret_cast<msComboInEq *>(data);
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
  if (!grad.empty()) {
    for (size_t i = 0; i < grad.size(); i++) {
      grad[i] = 0.0;
    }
    for (int l = 0; l < nT; l++) {
      m = m - 1;
      double iTop = degree[l] + iOffset;
      double iBottom = iOffset;
      for (int k = 0; k < nObs[l]; k++) {
        double sum = x[iTop];
        for (int j = iTop - 1; j >= iBottom; j--) {
          sum = sum * doses[m][k] + x[j];
        }
        if (sum < 0) sum = 0.0;
        double P = 1.0 - exp(-1.0 * sum);
        resid = (Y[m][k] * dslog(P) - (n_group[m][k] - Y[m][k]) * dslog(1.0 - P)) * (P - 1.0);
        int iIndex = iTop;
        for (int j = degree[l]; j >= 0; j--) {
          grad[iIndex] = grad[iIndex] + resid * (pow(doses[m][k], j));
          iIndex = iIndex - 1;
        }
      }
      iOffset = iTop + 1;
    }
  }

  double sum2 = 0;
  iOffset = 1;
  m = nT;
  for (int l = 0; l < nT; l++) {
    m = m - 1;
    double iTop = degree[l] + iOffset;
    double iBottom = iOffset;
    for (int k = 0; k < nObs[l]; k++) {
      double sum = x[iTop];
      for (int j = iTop - 1; j >= iBottom; j--) {
        sum = sum * doses[m][k] + x[j];
      }
      if (sum < 0) sum = 0.0;
      double P = 1.0 - exp(-1 * sum);
      sum2 = sum2 + Y[m][k] * slog(P) + (n_group[m][k] - Y[m][k]) * slog(1.0 - P);
    }
    iOffset = iTop + 1;
  }
  return target - sum2;
}

double slog(double X) {
  double coefs[4] = {6.7165863851209542e50, -2.0154759155362862e+35, 2.0169759155362859e+19, -710};
  if (X >= 1e-16) {
    return log(X);
  } else {
    double v = 0.0;
    for (int i = 0; i < 4; i++) {
      v = X * v + coefs[i];
    }
    return v;
  }
}

double dslog(double P) {
  if (P >= 1e-10) {
    return 1.0 / P;
  } else {
    return 2.0e10 - 1.0e20 * P;
  }
}

void BMDS_ENTRY_API __stdcall runBMDSDichoAnalysis(
    struct dichotomous_analysis *anal, struct dichotomous_model_result *res,
    struct dichotomous_GOF *gof, struct BMDS_results *bmdsRes, struct dicho_AOD *bmdsAOD
) {
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
  double *gofExpected = (double *)malloc(anal->n * sizeof(double));
  double *gofResidual = (double *)malloc(anal->n * sizeof(double));
  double gofTestStat = BMDS_MISSING;
  double gofPVal = BMDS_MISSING;
  double gofDF = BMDS_MISSING;
  double *ebUpper = (double *)malloc(anal->n * sizeof(double));
  double *ebLower = (double *)malloc(anal->n * sizeof(double));
  gofRes.n = anal->n;
  gofRes.expected = gofExpected;
  gofRes.residual = gofResidual;
  gofRes.test_statistic = gofTestStat;

  gofRes.p_value = gofPVal;
  gofRes.df = gofDF;

  compute_dichotomous_pearson_GOF(&gofData, &gofRes);

  gof->test_statistic = gofRes.test_statistic;

  // these will be updated later with bounded parm info
  gof->p_value = gofRes.p_value;
  gof->df = gofRes.df;

  gof->n = gofRes.n;
  for (int i = 0; i < gofRes.n; i++) {
    gof->expected.push_back(gofRes.expected[i]);
    gof->residual.push_back(gofRes.residual[i]);
  }

  // do error bar calcs
  double pHat;
  double z;
  double eb1;
  double eb2;
  double ebDenom;
  double gofAlpha = 0.05;                          // Alpha value for 95% confidence limit
  z = gsl_cdf_ugaussian_Pinv(1.0 - gofAlpha / 2);  // Z score
  for (int i = 0; i < gof->n; i++) {
    pHat = anal->Y[i] / anal->n_group[i];  // observed probability
    eb1 = (2 * anal->Y[i] + z * z);
    eb2 =
        z *
        sqrt(z * z - (2 + 1 / anal->n_group[i]) + 4 * pHat * ((anal->n_group[i] - anal->Y[i]) + 1));
    ebDenom = 2.0 * (anal->n_group[i] + z * z);
    gof->ebLower.push_back((eb1 - 1 - eb2) / ebDenom);
    gof->ebUpper.push_back((eb1 + 1 + eb2) / ebDenom);
  }

  // calculate model chi^2 value
  bmdsRes->chisq = 0.0;
  for (int i = 0; i < gofRes.n; i++) {
    bmdsRes->chisq += gofRes.residual[i] * gofRes.residual[i];
  }

  // calculate bayesian BIC_equiv
  Eigen::MatrixXd cov(res->nparms, res->nparms);
  int row = 0;
  int col = 0;
  for (int i = 0; i < res->nparms * res->nparms; i++) {
    col = i / res->nparms;
    row = i - col * res->nparms;
    cov(row, col) = res->cov[i];
  }

  bmdsRes->BIC_equiv =
      res->nparms / 2.0 * log(2.0 * M_PI) + res->max + 0.5 * log(max(0.0, cov.determinant()));
  bmdsRes->BIC_equiv = -1 * bmdsRes->BIC_equiv;

  // calculate dichtomous analysis of deviance
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

  // incorporate affect of bounded parameters
  int bounded = 0;
  for (int i = 0; i < anal->parms; i++) {
    if (bmdsRes->bounded[i]) {
      bounded++;
    }
  }

  bmdsAOD->nFit = anal->parms - bounded;     // number of estimated parameter
  bmdsAOD->dfFit = anal->n - bmdsAOD->nFit;  // nObs - nEstParms

  if (bmdsAOD->devFit < 0 || bmdsAOD->dfFit < 0) {
    bmdsAOD->pvFit = BMDS_MISSING;
  } else {
    bmdsAOD->pvFit = 1.0 - gsl_cdf_chisq_P(bmdsAOD->devFit, bmdsAOD->dfFit);
  }

  // update df for frequentist models only
  if (anal->prior[1] == 0) {
    // frequentist
    gofRes.df = bmdsAOD->dfFit;
  }
  gof->df = gofRes.df;

  if (gof->df > 0.0) {
    gofRes.p_value = 1.0 - gsl_cdf_chisq_P(bmdsRes->chisq, gof->df);
  } else {
    gofRes.p_value = 1.0;
  }

  if (gof->df <= 0.0) {
    gof->p_value = BMDS_MISSING;
  } else {
    gof->p_value = gofRes.p_value;
  }

  for (int i = 0; i < anal->parms; i++) {
    // std err is sqrt of covariance diagonals unless parameter hit a bound, then report NA
    bmdsRes->stdErr.push_back(BMDS_MISSING);
    bmdsRes->lowerConf.push_back(BMDS_MISSING);
    bmdsRes->upperConf.push_back(BMDS_MISSING);
  }

  calcParmCIs_dicho(res, bmdsRes);

  // compare Matt's BMD to BMD from CDF
  if (abs(bmdsRes->BMD - res->bmd) > BMDS_EPS) {
    std::cout << "Warning: CDF BMD differs from model BMD" << std::endl;
  }
  // Use Matt's BMD by default
  bmdsRes->BMD = res->bmd;

  clean_dicho_results(res, gof, bmdsRes, bmdsAOD);

  bool goodCDF = false;
  for (int i = 0; i < res->dist_numE; i++) {
    if (res->bmd_dist[i] != BMDS_MISSING) {
      goodCDF = true;
    }
  }

  if (bmdsRes->BMD != BMDS_MISSING && !std::isinf(bmdsRes->BMD) && std::isfinite(bmdsRes->BMD) &&
      goodCDF) {
    bmdsRes->validResult = true;
  }
}

void BMDS_ENTRY_API __stdcall runBMDSDichoMA(
    struct dichotomousMA_analysis *MA, struct dichotomous_analysis *DA,
    struct dichotomousMA_result *res, struct BMDSMA_results *bmdsRes
) {
  bmdsRes->BMD_MA = BMDS_MISSING;
  bmdsRes->BMDL_MA = BMDS_MISSING;
  bmdsRes->BMDU_MA = BMDS_MISSING;

  estimate_ma_laplace_dicho(MA, DA, res);

  collect_dichoMA_bmd_values(MA, res, bmdsRes, DA->alpha);
  for (int i = 0; i < MA->nmodels; i++) {
    rescale_dichoParms(MA->models[i], res->models[i]->parms);
  }

  // calc error bars for graphs
  double pHat;
  double z;
  double eb1;
  double eb2;
  double ebDenom;
  double gofAlpha = 0.05;                          // Alpha value for 95% confidence limit
  z = gsl_cdf_ugaussian_Pinv(1.0 - gofAlpha / 2);  // Z score
  for (int i = 0; i < DA->n; i++) {
    pHat = DA->Y[i] / DA->n_group[i];  // observed probability
    eb1 = (2 * DA->Y[i] + z * z);
    eb2 = z * sqrt(z * z - (2 + 1 / DA->n_group[i]) + 4 * pHat * ((DA->n_group[i] - DA->Y[i]) + 1));
    ebDenom = 2.0 * (DA->n_group[i] + z * z);
    bmdsRes->ebLower[i] = (eb1 - 1 - eb2) / ebDenom;
    bmdsRes->ebUpper[i] = (eb1 + 1 + eb2) / ebDenom;
  }

  clean_dicho_MA_results(res, bmdsRes);
}

// transform parameters or priors which were computed using a logistic dist.
void rescale_dichoParms(int model, double *parms) {
  // rescale background parameters for all models and v parameter for DHill model
  switch (model) {
    case dich_model::d_multistage:
    case dich_model::d_weibull:
    case dich_model::d_gamma:
    case dich_model::d_loglogistic:
    case dich_model::d_qlinear:
    case dich_model::d_logprobit:
      // rescale background parameter
      parms[0] = 1.0 / (1.0 + exp(-1.0 * parms[0]));
      break;
    case dich_model::d_hill:
      // rescale background and v parameter
      parms[0] = 1.0 / (1.0 + exp(-1.0 * parms[0]));
      parms[1] = 1.0 / (1.0 + exp(-1.0 * parms[1]));
      break;
    default:
      break;
  }
}

void calcParmCIs_dicho(struct dichotomous_model_result *res, struct BMDS_results *bmdsRes) {
  double *adj = new double[res->nparms];  // adjustments for transformed parameters
  int diagInd = 0;                        // 1d index of 2d diagonal location
  double tmp = 0;                         // used for tmp calculations

  for (int i = 0; i < res->nparms; i++) {
    diagInd = i * (1 + res->nparms);  // gives index of diagonal for 1d flattened array
    adj[i] = 1.0;
    if (!bmdsRes->bounded[i]) {
      bmdsRes->stdErr[i] = sqrt(res->cov[diagInd]);
    }
  }

  // now check if parameters were transformed and modify adj values accordingly
  // if transformed, then multiply adj * derivative of transform (at parm value)
  switch (res->model) {
    // transform parameters which were computed using a logistic dist.
    // void rescale_dichoParms(int model, struct dichotomous_model_result *res){
    case dich_model::d_multistage:
    case dich_model::d_weibull:
    case dich_model::d_gamma:
    case dich_model::d_loglogistic:
    case dich_model::d_qlinear:
    case dich_model::d_logprobit:
      // background parameter
      tmp = exp(-1 * res->parms[0]);
      tmp = tmp / ((1.0 + tmp) * (1.0 + tmp));
      adj[0] = tmp * tmp;
      break;
    case dich_model::d_hill:
      // background and v parameter
      tmp = exp(-1 * res->parms[0]);
      tmp = tmp / ((1.0 + tmp) * (1.0 + tmp));
      adj[0] = tmp * tmp;
      tmp = exp(-1 * res->parms[1]);
      tmp = tmp / ((1.0 + tmp) * (1.0 + tmp));
      adj[1] = tmp * tmp;
      break;
    default:
      break;
  }

  // now calculate upper and lower confidence intervals

  for (int i = 0; i < res->nparms; i++) {
    if (bmdsRes->bounded[i]) {
      bmdsRes->lowerConf[i] = BMDS_MISSING;
      bmdsRes->upperConf[i] = BMDS_MISSING;
    } else {
      bmdsRes->stdErr[i] = bmdsRes->stdErr[i] * adj[i];
      bmdsRes->lowerConf[i] = res->parms[i] - bmdsRes->stdErr[i] * BMDS_QNORM;
      bmdsRes->upperConf[i] = res->parms[i] + bmdsRes->stdErr[i] * BMDS_QNORM;
    }
  }

  delete[] adj;
}

void BMDS_ENTRY_API __stdcall runBMDSContAnalysis(
    struct continuous_analysis *anal, struct continuous_model_result *res,
    struct BMDS_results *bmdsRes, struct continuous_AOD *aod, struct continuous_GOF *gof,
    bool *detectAdvDir, bool *restricted
) {
  bmdsRes->BMD = BMDS_MISSING;
  bmdsRes->BMDL = BMDS_MISSING;
  bmdsRes->BMDU = BMDS_MISSING;
  bmdsRes->validResult = false;
  anal->transform_dose = false;
  if (anal->model != cont_model::exp_3 && anal->model != cont_model::exp_5) {
    if (anal->disttype == distribution::log_normal) {
      return;
    }
  }
  if (*detectAdvDir) {
    determineAdvDir(anal);
    int ind;
    if (*restricted) {
      switch (anal->model) {
        case cont_model::exp_3:
        case cont_model::exp_5:
          if (anal->prior[0] == 0) {  // checks if frequentist model
            if (anal->isIncreasing) {
              if (anal->disttype == distribution::normal_ncv) {
                anal->prior[20] = 0.0;  // c min
              } else {
                anal->prior[17] = 0.0;  // c min
              }
            } else {
              if (anal->disttype == distribution::normal_ncv) {
                anal->prior[26] = 0.0;  // c max
              } else {
                anal->prior[22] = 0.0;  // c max
              }
            }
          }
          break;
        case cont_model::polynomial:
          if (anal->prior[0] == 0) {  // checks if frequentist model
            int numRows = 2 + anal->degree;
            int ind;  // index in prior array
            if (anal->disttype == distribution::normal_ncv) {
              numRows++;
            }
            // set ind to target min or max for 1st beta parameter then adjust as needed
            // min for increasing, max for decreasing
            if (anal->isIncreasing) {
              ind = 3 * numRows + 1;
            } else {
              ind = 4 * numRows + 1;
            }
            // beta terms
            for (int i = ind; i < ind + anal->degree; i++) {
              anal->prior[i] = 0.0;  // beta min or max
            }
          }
          break;
      }  // end switch
    }  // end if restricted
  }  // end if detectAdvDir

  // need to handle issue with decreasing response datasets for relative deviation
  // easiest workaround is to modify BMRF before sending to estimate_sm_laplace
  if (anal->BMD_type == 3) {  // rel dev
    if (!anal->isIncreasing) {
      anal->BMR = 1.0 - anal->BMR;
    }
  }

  estimate_sm_laplace_cont(anal, res);

  // if not suff_stat, then convert
  struct continuous_analysis GOFanal;
  // arrays are needed for conversion to suff_stat
  double *doses = (double *)malloc(anal->n * sizeof(double));
  double *means = (double *)malloc(anal->n * sizeof(double));
  double *n_group = (double *)malloc(anal->n * sizeof(double));
  double *sd = (double *)malloc(anal->n * sizeof(double));
  for (int i = 0; i < anal->n; i++) {
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

  if (anal->suff_stat) {
    GOFanal = *anal;
  } else {
    // copy analysis and convert to suff_stat
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
    GOFanal.disttype = distribution::normal;  // needed for all distrubutions to avoid error in ln
                                              // conversion to suff_stats
    bmdsConvertSStat(anal, &GOFanal, true);
  }

  continuous_expected_result GOFres;
  GOFres.n = GOFanal.n;
  GOFres.expected = new double[GOFanal.n];
  GOFres.sd = new double[GOFanal.n];

  continuous_expectation(&GOFanal, res, &GOFres);

  for (int i = 0; i < GOFanal.n; i++) {
    gof->dose.push_back(GOFanal.doses[i]);
    gof->size.push_back(GOFanal.n_group[i]);
    gof->estMean.push_back(GOFres.expected[i]);
    gof->obsMean.push_back(GOFanal.Y[i]);
    gof->estSD.push_back(GOFres.sd[i]);
    gof->obsSD.push_back(GOFanal.sd[i]);
    gof->res.push_back(sqrt(gof->size[i]) * (gof->obsMean[i] - gof->estMean[i]) / gof->estSD[i]);
  }
  gof->n = GOFanal.n;
  if (anal->disttype == distribution::log_normal) {
    for (int i = 0; i < GOFanal.n; i++) {
      gof->calcMean.push_back(
          exp(log(GOFanal.Y[i]) - log(1 + pow(GOFanal.sd[i] / GOFanal.Y[i], 2.0)) / 2)
      );
      gof->calcSD.push_back(exp(sqrt(log(1.0 + pow(GOFanal.sd[i] / GOFanal.Y[i], 2.0)))));
    }
    for (int i = 0; i < GOFanal.n; i++) {
      gof->estMean.push_back(exp(gof->estMean[i] + pow(exp(res->parms[res->nparms - 1]), 2) / 2));
      gof->res.push_back(sqrt(gof->size[i]) * (gof->obsMean[i] - gof->estMean[i]) / gof->estSD[i]);
    }
  } else {
    for (int i = 0; i < GOFanal.n; i++) {
      gof->calcMean.push_back(GOFanal.Y[i]);
      gof->calcSD.push_back(GOFanal.sd[i]);
    }
  }

  double ebUpper, ebLower;
  for (int i = 0; i < GOFanal.n; i++) {
    ebLower = gof->calcMean[i] +
              gsl_cdf_tdist_Pinv(0.025, gof->n - 1) * (gof->obsSD[i] / sqrt(gof->size[i]));
    ebUpper = gof->calcMean[i] +
              gsl_cdf_tdist_Pinv(0.975, gof->n - 1) * (gof->obsSD[i] / sqrt(gof->size[i]));
    gof->ebLower.push_back(ebLower);
    gof->ebUpper.push_back(ebUpper);
  }
  // calculate bayesian BIC_equiv
  Eigen::MatrixXd cov(res->nparms, res->nparms);
  int row = 0;
  int col = 0;
  for (int i = 0; i < res->nparms * res->nparms; i++) {
    col = i / res->nparms;
    row = i - col * res->nparms;
    cov(row, col) = res->cov[i];
  }

  bmdsRes->BIC_equiv =
      res->nparms / 2.0 * log(2.0 * M_PI) + res->max + 0.5 * log(max(0.0, cov.determinant()));
  bmdsRes->BIC_equiv = -1 * bmdsRes->BIC_equiv;
  collect_cont_bmd_values(anal, res, bmdsRes);

  aod->LL.resize(5);
  aod->nParms.resize(5);
  aod->AIC.resize(5);
  aod->TOI.llRatio.resize(4);
  aod->TOI.DF.resize(4);
  aod->TOI.pVal.resize(4);
  calc_contAOD(anal, &GOFanal, res, bmdsRes, aod);

  rescale_contParms(anal, res->parms);

  for (int i = 0; i < anal->parms; i++) {
    bmdsRes->stdErr.push_back(BMDS_MISSING);
    bmdsRes->lowerConf.push_back(BMDS_MISSING);
    bmdsRes->upperConf.push_back(BMDS_MISSING);
  }

  calcParmCIs_cont(res, bmdsRes);

  // compare Matt's BMD to BMD from CDF
  if (abs(bmdsRes->BMD - res->bmd) > BMDS_EPS) {
    std::cout << "Warning: CDF BMD differs from model BMD" << std::endl;
  }
  // Use Matt's BMD by default
  bmdsRes->BMD = res->bmd;

  clean_cont_results(res, bmdsRes, aod, gof);

  bool goodCDF = false;
  for (int i = 0; i < res->dist_numE; i++) {
    if (res->bmd_dist[i] != BMDS_MISSING) {
      goodCDF = true;
    }
  }

  if (bmdsRes->BMD != BMDS_MISSING && !std::isinf(bmdsRes->BMD) && std::isfinite(bmdsRes->BMD) &&
      goodCDF) {
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
  delete[] GOFres.expected;
  delete[] GOFres.sd;
}

void rescale_contParms(struct continuous_analysis *CA, double *parms) {
  // rescale alpha parameter for all continuous models
  // assumes alpha parameter is always last
  switch (CA->model) {
    case cont_model::hill:
    case cont_model::power:
    case cont_model::funl:
    case cont_model::polynomial:
      // rescale alpha
      parms[CA->parms - 1] = exp(parms[CA->parms - 1]);
      break;
    case cont_model::exp_5:
      parms[2] = exp(parms[2]);
      break;
  }
}

void calcParmCIs_cont(struct continuous_model_result *res, struct BMDS_results *bmdsRes) {
  double *adj = new double[res->nparms];  // adjustments for transformed parameters
  int diagInd = 0;                        // 1d index of 2d diagonal location
  double tmp = 0;                         // used for tmp calculations

  for (int i = 0; i < res->nparms; i++) {
    diagInd = i * (1 + res->nparms);  // gives index of diagonal for 1d flattened array
    adj[i] = 1.0;
    if (!bmdsRes->bounded[i]) {
      bmdsRes->stdErr[i] = sqrt(res->cov[diagInd]);
    }
  }

  // now check if parameters were transformed and modify adj values accordingly
  // if transformed, then multiply adj * derivative of transform (at parm value)
  switch (res->model) {
    case cont_model::hill:
    case cont_model::power:
    case cont_model::funl:
    case cont_model::polynomial:
      // rescale alpha  //at this point alpha has already been rescaled, so no need to use exp(parm)
      // since actual parm is ln(alpha), then gprime(theta) = exp(theta) = exp(ln(alpha)) = alpha
      tmp = res->parms[res->nparms - 1];
      adj[res->nparms - 1] = adj[res->nparms - 1] * tmp * tmp;
      break;
    case cont_model::exp_5:
      // rescale c
      // tmp = exp(res->parms[2]);
      // adj[2] = adj[2]*tmp*tmp;
      break;
    default:
      break;
  }

  // now calculate upper and lower confidence intervals
  for (int i = 0; i < res->nparms; i++) {
    if (bmdsRes->bounded[i]) {
      bmdsRes->lowerConf[i] = BMDS_MISSING;
      bmdsRes->upperConf[i] = BMDS_MISSING;
    } else {
      bmdsRes->stdErr[i] = bmdsRes->stdErr[i] * adj[i];
      bmdsRes->lowerConf[i] = res->parms[i] - bmdsRes->stdErr[i] * BMDS_QNORM;
      bmdsRes->upperConf[i] = res->parms[i] + bmdsRes->stdErr[i] * BMDS_QNORM;
    }
  }
}

void calc_dichoAOD(
    struct dichotomous_analysis *DA, struct dichotomous_model_result *res,
    struct BMDS_results *bmdsRes, struct dicho_AOD *bmdsAOD, struct dichotomous_aod *aod
) {
  deviance_dichotomous(DA, aod);

  bmdsAOD->fullLL = -1 * aod->A1;
  bmdsAOD->nFull = aod->N1;
  bmdsAOD->redLL = -1 * aod->A2;
  bmdsAOD->nRed = aod->N2;
  bmdsAOD->fittedLL = -1 * res->max;
  int bounded = 0;
  for (int i = 0; i < DA->parms; i++) {
    if (bmdsRes->bounded[i]) {
      bounded++;
    }
  }

  bmdsAOD->devFit = 2 * (bmdsAOD->fullLL - bmdsAOD->fittedLL);

  // these fitted model values will be calculated later after bounded parms have been determined
  bmdsAOD->nFit = BMDS_MISSING;
  bmdsAOD->dfFit = BMDS_MISSING;
  bmdsAOD->pvFit = BMDS_MISSING;

  double dev;
  double df;

  bmdsAOD->devRed = dev = 2 * (bmdsAOD->fullLL - bmdsAOD->redLL);
  bmdsAOD->dfRed = df = DA->n - 1;
  if (dev < 0 || df < 0) {
    bmdsAOD->pvRed = BMDS_MISSING;
  } else {
    bmdsAOD->pvRed = 1.0 - gsl_cdf_chisq_P(dev, df);
  }
}

void calc_contAOD(
    struct continuous_analysis *CA, struct continuous_analysis *GOFanal,
    struct continuous_model_result *res, struct BMDS_results *bmdsRes, struct continuous_AOD *aod
) {
  continuous_deviance CD;

  if (CA->disttype == distribution::log_normal) {
    estimate_log_normal_aod(CA, &CD);
  } else {
    estimate_normal_aod(CA, &CD);
  }
  // fill aod with results and calculate tests of interest
  aod->LL[0] = -1 * CD.A1;
  aod->nParms[0] = CD.N1;
  aod->LL[1] = -1 * CD.A2;
  aod->nParms[1] = CD.N2;
  // A3 is strictly NCV test.  If CV use results from A1
  if (CA->disttype == distribution::normal_ncv) {
    aod->LL[2] = -1 * CD.A3;
    aod->nParms[2] = CD.N3;
  } else {
    aod->LL[2] = -1 * CD.A1;
    aod->nParms[2] = CD.N1;
  }
  aod->LL[3] = -1 * res->max;
  int bounded = 0;
  for (int i = 0; i < CA->parms; i++) {
    if (bmdsRes->bounded[i]) {
      bounded++;
    }
  }
  aod->nParms[3] = res->nparms - bounded;
  aod->LL[4] = -1 * CD.R;
  aod->nParms[4] = CD.NR;

  // add tests of interest
  // TODO:  need to check for bounded parms in A1,A2,A3,R
  double sumN = 0;
  for (int i = 0; i < 5; i++) {
    aod->AIC[i] = 2 * (-1 * aod->LL[i] + aod->nParms[i]);
  }

  if (CA->suff_stat) {
    for (int i = 0; i < CA->n; i++) {
      sumN += CA->n_group[i];
    }
  } else {
    sumN = CA->n;
  }
  aod->addConst = -1 * sumN * log(2 * M_PI) / 2;
  if (CA->disttype == distribution::log_normal) {
    double tmp = 0;
    if (CA->suff_stat) {
      for (int i = 0; i < CA->n; i++) {
        tmp += (log(CA->Y[i]) - log(1 + pow((CA->sd[i] / CA->Y[i]), 2)) / 2) * CA->n_group[i];
      }
    } else {
      GOFanal->disttype =
          distribution::log_normal;          // previous GOFanal is always normal distribution
      bmdsConvertSStat(CA, GOFanal, false);  // recalculate using with lognormal transformation
      double divisorTerm = GOFanal->Y[0];
      // rescale to match BMDS 3.x
      for (int i = 0; i < GOFanal->n; i++) {
        GOFanal->Y[i] /= divisorTerm;
        GOFanal->sd[i] /= divisorTerm;
      }
      // convert from logscale
      double tmpMean, tmpSD;
      for (int i = 0; i < GOFanal->n; i++) {
        tmpMean = GOFanal->Y[i];
        tmpSD = GOFanal->sd[i];
        GOFanal->Y[i] = exp(tmpMean + pow(tmpSD, 2) / 2);
        GOFanal->sd[i] = GOFanal->Y[i] * sqrt(exp(pow(tmpSD, 2)) - 1);

        // reverse response scaling
        GOFanal->Y[i] *= divisorTerm;
        GOFanal->sd[i] *= divisorTerm;

        // convert to log scale
        tmpMean = GOFanal->Y[i];
        tmpSD = GOFanal->sd[i];
        GOFanal->Y[i] = log(tmpMean) - log(1 + pow(tmpSD / tmpMean, 2)) / 2;
        GOFanal->sd[i] = sqrt(log(1 + pow(tmpSD / tmpMean, 2)));

        tmp += GOFanal->Y[i] * GOFanal->n_group[i];
      }
    }
    aod->addConst -= tmp;
  }

  calcTestsOfInterest(aod);
}

void calcTestsOfInterest(struct continuous_AOD *aod) {
  double dev;
  int df;

  // Test #1 - A2 vs Reduced - does mean and/or variance differ across dose groups
  dev = (aod->LL[1] - aod->LL[4]);
  // handle underflow/negative zero
  if (dev < BMDS_EPS) {
    dev = 0.0;
  }
  aod->TOI.llRatio[0] = dev = 2 * dev;

  aod->TOI.DF[0] = df = aod->nParms[1] - aod->nParms[4];
  aod->TOI.pVal[0] =
      (std::isnan(dev) || dev < 0.0 || df < 0.0) ? -1 : 1.0 - gsl_cdf_chisq_P(dev, df);

  // Test #2 - A1 vs A2 - homogeneity of variance across dose groups
  dev = (aod->LL[1] - aod->LL[0]);
  // handle underflow/negative zero
  if (dev < BMDS_EPS) {
    dev = 0.0;
  }
  aod->TOI.llRatio[1] = dev = 2 * dev;
  aod->TOI.DF[1] = df = aod->nParms[1] - aod->nParms[0];
  aod->TOI.pVal[1] =
      (std::isnan(dev) || dev < 0.0 || df < 0.0) ? -1 : 1.0 - gsl_cdf_chisq_P(dev, df);

  // Test #3 - A2 vs A3 - Does the model describe variances adequately
  dev = (aod->LL[1] - aod->LL[2]);
  // handle underflow/negative zero
  if (dev < BMDS_EPS) {
    dev = 0.0;
  }
  aod->TOI.llRatio[2] = dev = 2 * dev;
  aod->TOI.DF[2] = df = aod->nParms[1] - aod->nParms[2];
  aod->TOI.pVal[2] =
      (std::isnan(dev) || dev < 0.0 || df < 0.0) ? -1 : 1.0 - gsl_cdf_chisq_P(dev, df);

  // Test #4 - A3 vs Fitted - does the fitted model describe the obs data adequately
  dev = (aod->LL[2] - aod->LL[3]);
  // handle underflow/negative zero
  if (dev < BMDS_EPS) {
    dev = 0.0;
  }
  aod->TOI.llRatio[3] = dev = 2 * dev;
  aod->TOI.DF[3] = df = aod->nParms[2] - aod->nParms[3];
  aod->TOI.pVal[3] =
      (std::isnan(dev) || dev < 0.0 || df < 0.0) ? -1 : 1.0 - gsl_cdf_chisq_P(dev, df);

  // DOF check for test 4
  if (aod->TOI.DF[3] <= 0) {
    aod->TOI.pVal[3] = BMDS_MISSING;
  }
}

void determineAdvDir(struct continuous_analysis *CA) {
  int n_rows;

  struct continuous_analysis CAnew;
  // allocate memory for individual size.  Will not use entire array for summary data.
  double *tmpD = (double *)malloc(CA->n * sizeof(double));
  double *tmpY = (double *)malloc(CA->n * sizeof(double));
  double *tmpN = (double *)malloc(CA->n * sizeof(double));
  double *tmpSD = (double *)malloc(CA->n * sizeof(double));

  CAnew.doses = tmpD;
  CAnew.Y = tmpY;
  CAnew.n_group = tmpN;
  CAnew.sd = tmpSD;

  if (!CA->suff_stat) {
    bmdsConvertSStat(CA, &CAnew, true);
    n_rows = CAnew.n;
  } else {
    // already suff stat
    n_rows = CA->n;
  }
  Eigen::MatrixXd X(n_rows, 2);
  Eigen::MatrixXd W(n_rows, n_rows);
  Eigen::VectorXd Y(n_rows);
  Eigen::VectorXd beta(2);

  if (!CA->suff_stat) {
    for (int i = 0; i < n_rows; i++) {
      X(i, 0) = 1;
      X(i, 1) = CAnew.doses[i];

      W(i, i) = CAnew.n_group[i];
      Y(i) = CAnew.Y[i];
    }

  } else {
    for (int i = 0; i < n_rows; i++) {
      X(i, 0) = 1;
      X(i, 1) = CA->doses[i];

      W(i, i) = CA->n_group[i];
      Y(i) = CA->Y[i];
    }
  }

  beta = (X.transpose() * W * X).colPivHouseholderQr().solve(X.transpose() * W * Y);

  if (beta(1) > 0) {
    CA->isIncreasing = true;
  } else {
    CA->isIncreasing = false;
  }

  free(tmpD);
  free(tmpY);
  free(tmpN);
  free(tmpSD);
}

void bmdsConvertSStat(
    struct continuous_analysis *CA, struct continuous_analysis *CAss, bool clean
) {
  // standardize the data
  int n_rows = CA->n;
  int n_cols = CA->suff_stat ? 3 : 1;

  if (!CA->suff_stat) {
    Eigen::MatrixXd Yin(n_rows, n_cols);
    Eigen::MatrixXd Xin(n_rows, 1);
    Eigen::MatrixXd SSTAT, SSTAT_LN, X, Y;
    // copy the original data
    for (int i = 0; i < n_rows; i++) {
      Yin(i, 0) = CA->Y[i];
      Xin(i, 0) = CA->doses[i];
    }
    bool canConvert = convertSStat(Yin, Xin, &SSTAT, &SSTAT_LN, &X);
    bool isLognormal = CAss->disttype == distribution::log_normal;
    if (canConvert) {
      if (clean) {
        if (isLognormal) {
          Y = cleanSuffStat(SSTAT_LN, X, true, false);
        } else {
          Y = cleanSuffStat(SSTAT, X, false, false);
        }
      } else {
        if (isLognormal) {
          Y = SSTAT_LN;
        } else {
          Y = SSTAT;
        }
      }
      n_rows = X.rows();

      for (int i = 0; i < n_rows; i++) {
        CAss->doses[i] = X(i);
        CAss->Y[i] = Y.col(0)[i];
        CAss->n_group[i] = Y.col(1)[i];
        CAss->sd[i] = Y.col(2)[i];
      }
      CAss->n = n_rows;
    }
  } else {
    // already suff_stat, so just copy data to arrays
    for (int i = 0; i < n_rows; i++) {
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

// replace any double values containing NaN or infinity with BMDS_MISSING
void clean_dicho_results(
    struct dichotomous_model_result *res, struct dichotomous_GOF *gof, struct BMDS_results *bmdsRes,
    struct dicho_AOD *aod
) {
  // dichotomous_model_result
  for (int i = 0; i < res->nparms; i++) {
    cleanDouble(&res->parms[i]);
  }

  for (int i = 0; i < res->nparms * res->nparms; i++) {
    cleanDouble(&res->cov[i]);
  }

  cleanDouble(&res->max);
  cleanDouble(&res->model_df);
  cleanDouble(&res->total_df);
  cleanDouble(&res->bmd);

  for (int i = 0; i < res->dist_numE * 2; i++) {
    cleanDouble(&res->bmd_dist[i]);
  }

  // dichotomous_GOF
  for (int i = 0; i < gof->n; i++) {
    cleanDouble(&gof->expected[i]);
    cleanDouble(&gof->residual[i]);
    cleanDouble(&gof->ebLower[i]);
    cleanDouble(&gof->ebUpper[i]);
  }
  cleanDouble(&gof->test_statistic);
  cleanDouble(&gof->p_value);
  cleanDouble(&gof->df);

  // BMDS_results
  cleanDouble(&bmdsRes->BMD);
  cleanDouble(&bmdsRes->BMDL);
  cleanDouble(&bmdsRes->BMDU);
  cleanDouble(&bmdsRes->AIC);
  cleanDouble(&bmdsRes->BIC_equiv);
  cleanDouble(&bmdsRes->chisq);

  for (int i = 0; i < res->nparms; i++) {
    cleanDouble(&bmdsRes->stdErr[i]);
    cleanDouble(&bmdsRes->lowerConf[i]);
    cleanDouble(&bmdsRes->upperConf[i]);
  }

  // dicho_AOD
  cleanDouble(&aod->fullLL);
  cleanDouble(&aod->redLL);
  cleanDouble(&aod->fittedLL);
  cleanDouble(&aod->devFit);
  cleanDouble(&aod->devRed);
  cleanDouble(&aod->pvFit);
  cleanDouble(&aod->pvRed);
}

void clean_cont_results(
    struct continuous_model_result *res, struct BMDS_results *bmdsRes, struct continuous_AOD *aod,
    struct continuous_GOF *gof
) {
  // continuous_model_result
  for (int i = 0; i < res->nparms; i++) {
    cleanDouble(&res->parms[i]);
  }
  for (int i = 0; i < res->nparms * res->nparms; i++) {
    cleanDouble(&res->cov[i]);
  }

  cleanDouble(&res->max);
  cleanDouble(&res->model_df);
  cleanDouble(&res->total_df);
  cleanDouble(&res->bmd);
  for (int i = 0; i < res->dist_numE * 2; i++) {
    cleanDouble(&res->bmd_dist[i]);
  }

  // BMDS_results
  cleanDouble(&bmdsRes->BMD);
  cleanDouble(&bmdsRes->BMDL);
  cleanDouble(&bmdsRes->BMDU);
  cleanDouble(&bmdsRes->AIC);
  cleanDouble(&bmdsRes->BIC_equiv);
  cleanDouble(&bmdsRes->chisq);
  for (int i = 0; i < res->nparms; i++) {
    cleanDouble(&bmdsRes->stdErr[i]);
    cleanDouble(&bmdsRes->lowerConf[i]);
    cleanDouble(&bmdsRes->upperConf[i]);
  }

  // continuous_AOD and testsOfInterest
  for (int i = 0; i < 5; i++) {
    cleanDouble(&aod->LL[i]);
    cleanDouble(&aod->AIC[i]);
    for (int j = 0; j < 4; j++) {
      cleanDouble(&aod->TOI.llRatio[j]);
      cleanDouble(&aod->TOI.DF[j]);
      cleanDouble(&aod->TOI.pVal[j]);
    }
  }
  cleanDouble(&aod->addConst);

  // continuous_GOF
  for (int i = 0; i < gof->n; i++) {
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
  // dichotomousMA_result
  for (int i = 0; i < res->nmodels; i++) {
    cleanDouble(&res->post_probs[i]);
    for (int j = 0; j < res->models[i]->nparms; j++) {
      cleanDouble(&res->models[i]->parms[j]);
    }
    for (int j = 0; j < res->models[i]->nparms * res->models[i]->nparms; j++) {
      cleanDouble(&res->models[i]->cov[j]);
    }
    cleanDouble(&res->models[i]->max);
    cleanDouble(&res->models[i]->model_df);
    cleanDouble(&res->models[i]->total_df);
    cleanDouble(&res->models[i]->bmd);
    for (int j = 0; j < res->models[i]->dist_numE * 2; j++) {
      cleanDouble(&res->models[i]->bmd_dist[j]);
    }
  }
  for (int i = 0; i < res->dist_numE * 2; i++) {
    cleanDouble(&res->bmd_dist[i]);
  }

  // BMDSMA_results
  cleanDouble(&bmdsRes->BMD_MA);
  cleanDouble(&bmdsRes->BMDL_MA);
  cleanDouble(&bmdsRes->BMDU_MA);
  for (int i = 0; i < res->nmodels; i++) {
    cleanDouble(&bmdsRes->BMD[i]);
    cleanDouble(&bmdsRes->BMDL[i]);
    cleanDouble(&bmdsRes->BMDU[i]);
  }
}

// replace any double values containing NaN or infinity with BMDS_MISSING
void cleanDouble(double *val) {
  if (!isfinite(*val)) {
    *val = BMDS_MISSING;
  }
}

string BMDS_ENTRY_API __stdcall version() { return BMDS_VERSION; }

void convertToPythonDichoRes(
    struct dichotomous_model_result *res, struct python_dichotomous_model_result *pyRes
) {
  pyRes->model = res->model;
  pyRes->nparms = res->nparms;
  pyRes->parms.assign(res->parms, res->parms + res->nparms);
  pyRes->cov.assign(res->cov, res->cov + res->nparms * res->nparms);
  pyRes->max = res->max;
  pyRes->dist_numE = res->dist_numE;
  pyRes->model_df = res->model_df;
  pyRes->total_df = res->total_df;
  pyRes->bmd_dist.assign(res->bmd_dist, res->bmd_dist + res->dist_numE * 2);
  pyRes->bmd = res->bmd;
  pyRes->gof_p_value = res->gof_p_value;
  pyRes->gof_chi_sqr_statistic = res->gof_chi_sqr_statistic;
}

void convertFromPythonDichoAnalysis(
    struct dichotomous_analysis *anal, struct python_dichotomous_analysis *pyAnal
) {
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

  if (pyAnal->n == pyAnal->doses.size() && pyAnal->doses.size() == pyAnal->Y.size() &&
      pyAnal->doses.size() == pyAnal->n_group.size()) {
    for (int i = 0; i < pyAnal->n; i++) {
      anal->Y[i] = pyAnal->Y[i];
      anal->doses[i] = pyAnal->doses[i];
      anal->n_group[i] = pyAnal->n_group[i];
    }
  }
  if (pyAnal->prior.size() > 0) {
    for (int i = 0; i < pyAnal->prior.size(); i++) {
      anal->prior[i] = pyAnal->prior[i];
    }
  }
}

void convertFromPythonDichoRes(
    struct dichotomous_model_result *res, struct python_dichotomous_model_result *pyRes
) {
  res->model = pyRes->model;
  res->nparms = pyRes->nparms;
  res->max = pyRes->max;
  res->dist_numE = pyRes->dist_numE;
  res->model_df = pyRes->model_df;
  res->total_df = pyRes->total_df;
  res->bmd = pyRes->bmd;
  res->gof_p_value = pyRes->gof_p_value;
  res->gof_chi_sqr_statistic = pyRes->gof_chi_sqr_statistic;
  // arrays & vectors
  if (pyRes->parms.size() > 0) {
    for (int i = 0; i < pyRes->parms.size(); i++) {
      res->parms[i] = pyRes->parms[i];
    }
  }
  if (pyRes->cov.size() > 0) {
    for (int i = 0; i < pyRes->cov.size(); i++) {
      res->cov[i] = pyRes->cov[i];
    }
  }
  if (pyRes->bmd_dist.size() > 0) {
    for (int i = 0; i < pyRes->bmd_dist.size(); i++) {
      res->bmd_dist[i] = pyRes->bmd_dist[i];
    }
  }
}

void convertFromPythonDichoMAAnalysis(
    struct dichotomousMA_analysis *MA, struct python_dichotomousMA_analysis *pyMA
) {
  MA->nmodels = pyMA->nmodels;
  for (int i = 0; i < pyMA->nmodels; i++) {
    MA->priors[i] = pyMA->priors[i].data();
  }
  MA->nparms = pyMA->nparms.data();
  MA->actual_parms = pyMA->actual_parms.data();
  MA->prior_cols = pyMA->prior_cols.data();
  MA->models = pyMA->models.data();
  MA->modelPriors = pyMA->modelPriors.data();
}

void convertFromPythonDichoMARes(
    struct dichotomousMA_result *res, struct python_dichotomousMA_result *pyRes
) {
  res->nmodels = pyRes->nmodels;
  for (int i = 0; i < pyRes->nmodels; i++) {
    convertFromPythonDichoRes(res->models[i], &pyRes->models[i]);
  }
  res->dist_numE = pyRes->dist_numE;
  res->post_probs = pyRes->post_probs.data();
  res->bmd_dist = pyRes->bmd_dist.data();
}

void convertFromPythonContAnalysis(
    struct continuous_analysis *anal, struct python_continuous_analysis *pyAnal
) {
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
  if (pyAnal->suff_stat) {
    validated = pyAnal->n == pyAnal->doses.size() && pyAnal->doses.size() == pyAnal->Y.size() &&
                pyAnal->doses.size() == pyAnal->n_group.size();
  } else {
    validated = pyAnal->n == pyAnal->doses.size() && pyAnal->doses.size() == pyAnal->Y.size();
  }

  if (validated) {
    for (int i = 0; i < pyAnal->n; i++) {
      anal->Y[i] = pyAnal->Y[i];
      anal->doses[i] = pyAnal->doses[i];
    }
  }
  if (validated && pyAnal->suff_stat) {
    for (int i = 0; i < pyAnal->n; i++) {
      anal->n_group[i] = pyAnal->n_group[i];
      anal->sd[i] = pyAnal->sd[i];
    }
  }
  if (pyAnal->prior.size() > 0) {
    for (int i = 0; i < pyAnal->prior.size(); i++) {
      anal->prior[i] = pyAnal->prior[i];
    }
  }
}

void convertFromPythonContRes(
    struct continuous_model_result *res, struct python_continuous_model_result *pyRes
) {
  res->model = pyRes->model;
  res->dist = pyRes->dist;
  res->nparms = pyRes->nparms;
  res->max = pyRes->max;
  res->dist_numE = pyRes->dist_numE;
  res->model_df = pyRes->model_df;
  res->total_df = pyRes->total_df;
  res->bmd = pyRes->bmd;
  // arrays & vectors
  if (pyRes->parms.size() > 0) {
    for (int i = 0; i < pyRes->parms.size(); i++) {
      res->parms[i] = pyRes->parms[i];
    }
  }
  if (pyRes->cov.size() > 0) {
    for (int i = 0; i < pyRes->cov.size(); i++) {
      res->cov[i] = pyRes->cov[i];
    }
  }
  if (pyRes->bmd_dist.size() > 0) {
    for (int i = 0; i < pyRes->bmd_dist.size(); i++) {
      res->bmd_dist[i] = pyRes->bmd_dist[i];
    }
  }
}

void convertToPythonContRes(
    struct continuous_model_result *res, struct python_continuous_model_result *pyRes
) {
  pyRes->model = res->model;
  pyRes->dist = res->dist;
  pyRes->nparms = res->nparms;
  pyRes->parms.assign(res->parms, res->parms + res->nparms);
  pyRes->cov.assign(res->cov, res->cov + res->nparms * res->nparms);
  pyRes->max = res->max;
  pyRes->dist_numE = res->dist_numE;
  pyRes->model_df = res->model_df;
  pyRes->total_df = res->total_df;
  pyRes->bmd_dist.assign(res->bmd_dist, res->bmd_dist + res->dist_numE * 2);
  pyRes->bmd = res->bmd;
}

void BMDS_ENTRY_API __stdcall pythonBMDSDicho(
    struct python_dichotomous_analysis *pyAnal, struct python_dichotomous_model_result *pyRes
) {
  // 1st convert from python struct
  dichotomous_analysis anal;
  anal.Y = new double[pyAnal->n];
  anal.doses = new double[pyAnal->n];
  anal.n_group = new double[pyAnal->n];
  anal.prior = new double[pyAnal->parms * pyAnal->prior_cols];
  convertFromPythonDichoAnalysis(&anal, pyAnal);

  dichotomous_model_result res;
  res.parms = new double[pyRes->nparms];
  res.cov = new double[pyRes->nparms * pyRes->nparms];
  res.bmd_dist = new double[pyRes->dist_numE * 2];
  convertFromPythonDichoRes(&res, pyRes);

  runBMDSDichoAnalysis(&anal, &res, &pyRes->gof, &pyRes->bmdsRes, &pyRes->aod);

  convertToPythonDichoRes(&res, pyRes);
}

void BMDS_ENTRY_API __stdcall pythonBMDSDichoMA(
    struct python_dichotomousMA_analysis *pyMA, struct python_dichotomousMA_result *pyRes
) {
  // convert python_dichtomousMA_analysis to dichotomousMA_analysis
  dichotomousMA_analysis MA;
  MA.priors = new double *[pyMA->nmodels];
  MA.nparms = new int[pyMA->nmodels];
  MA.actual_parms = new int[pyMA->nmodels];
  MA.prior_cols = new int[pyMA->nmodels];
  MA.models = new int[pyMA->nmodels];
  MA.modelPriors = new double[pyMA->nmodels];
  convertFromPythonDichoMAAnalysis(&MA, pyMA);

  // convert python_dichotomous_analysis to dichotomous_analysis
  dichotomous_analysis DA;
  DA.Y = new double[pyMA->pyDA.n];
  DA.doses = new double[pyMA->pyDA.n];
  DA.n_group = new double[pyMA->pyDA.n];
  DA.prior = new double[pyMA->pyDA.parms * pyMA->pyDA.prior_cols];
  convertFromPythonDichoAnalysis(&DA, &pyMA->pyDA);

  //  convertFromPythonDichoMARes(&res, pyRes);

  struct dichotomous_model_result **res = new struct dichotomous_model_result *[pyRes->nmodels];
  for (int i = 0; i < pyRes->nmodels; i++) {
    res[i] = (struct dichotomous_model_result *)malloc(sizeof(struct dichotomous_model_result));
    res[i]->model = pyRes->models[i].model;
    res[i]->nparms = pyRes->models[i].nparms;
    res[i]->dist_numE = pyRes->models[i].dist_numE;
    res[i]->parms = (double *)malloc(sizeof(double) * pyRes->models[i].nparms);
    res[i]->cov =
        (double *)malloc(sizeof(double) * pyRes->models[i].nparms * pyRes->models[i].nparms);
    res[i]->bmd_dist = (double *)malloc(sizeof(double) * pyRes->models[i].dist_numE * 2);
  }

  struct dichotomousMA_result maRes;
  maRes.nmodels = pyRes->nmodels;
  maRes.models = res;
  maRes.dist_numE = pyRes->dist_numE;
  maRes.post_probs = new double[pyRes->nmodels];
  maRes.bmd_dist = new double[pyRes->dist_numE * 2];

  runBMDSDichoMA(&MA, &DA, &maRes, &pyRes->bmdsRes);
  // convert back to python objs
  // pyDA should not change
  pyRes->post_probs.resize(pyRes->nmodels);
  for (int i = 0; i < pyRes->nmodels; i++) {
    pyRes->post_probs[i] = maRes.post_probs[i];
    pyRes->models[i].max = maRes.models[i]->max;
    pyRes->models[i].model_df = maRes.models[i]->model_df;
    pyRes->models[i].total_df = maRes.models[i]->total_df;
    pyRes->models[i].bmd = maRes.models[i]->bmd;
    pyRes->models[i].gof_p_value = maRes.models[i]->gof_p_value;
    pyRes->models[i].gof_chi_sqr_statistic = maRes.models[i]->gof_chi_sqr_statistic;
    int nparms = pyRes->models[i].nparms;

    pyRes->models[i].parms.resize(nparms);
    for (int j = 0; j < nparms; j++) {
      pyRes->models[i].parms[j] = maRes.models[i]->parms[j];
    }
    pyRes->models[i].cov.resize(nparms * nparms);
    for (int j = 0; j < nparms * nparms; j++) {
      pyRes->models[i].cov[j] = maRes.models[i]->cov[j];
    }
    pyRes->models[i].bmd_dist.resize(pyRes->dist_numE * 2);
    for (int j = 0; j < pyRes->dist_numE * 2; j++) {
      pyRes->models[i].bmd_dist[j] = maRes.models[i]->bmd_dist[j];
    }
  }
  pyRes->bmd_dist.resize(pyRes->dist_numE * 2);
  for (int i = 0; i < pyRes->dist_numE * 2; i++) {
    pyRes->bmd_dist[i] = maRes.bmd_dist[i];
  }
}

void BMDS_ENTRY_API __stdcall pythonBMDSCont(
    struct python_continuous_analysis *pyAnal, struct python_continuous_model_result *pyRes
) {
  // convert pyAnal to anal
  continuous_analysis anal;
  anal.Y = new double[pyAnal->n];
  anal.doses = new double[pyAnal->n];
  anal.sd = new double[pyAnal->n];
  anal.n_group = new double[pyAnal->n];
  anal.prior = new double[pyAnal->parms * pyAnal->prior_cols];
  convertFromPythonContAnalysis(&anal, pyAnal);

  // convert pyRes to res
  continuous_model_result res;
  res.parms = new double[pyRes->nparms];
  res.cov = new double[pyRes->nparms * pyRes->nparms];
  res.bmd_dist = new double[pyRes->dist_numE * 2];
  convertFromPythonContRes(&res, pyRes);

  runBMDSContAnalysis(
      &anal, &res, &pyRes->bmdsRes, &pyRes->aod, &pyRes->gof, &pyAnal->detectAdvDir,
      &pyAnal->restricted
  );

  convertToPythonContRes(&res, pyRes);
}

void BMDS_ENTRY_API __stdcall pythonBMDSMultitumor(
    struct python_multitumor_analysis *pyAnal, struct python_multitumor_result *pyRes
) {
  // run each individual multistage model
  for (int i = 0; i < pyAnal->ndatasets; i++) {
    // only all individual multistage models if degree == 0
    // else only run the specified model
    if (pyAnal->degree[i] == 0) {
      for (int j = 0; j < pyAnal->nmodels[i]; j++) {
        pythonBMDSDicho(&pyAnal->models[i][j], &pyRes->models[i][j]);
      }
    } else {
      pythonBMDSDicho(&pyAnal->models[i][0], &pyRes->models[i][0]);
    }
  }

  // select models
  selectMultitumorModel(pyAnal, pyRes);

  // create new pyAnal and pyRes only containing selected models
  // Note: ndatasets may differ from previous structs (pyAnal & pyRes) depending on whether any
  // datasets were rejected
  struct python_multitumor_analysis anal;
  struct python_multitumor_result res;

  anal.n = pyAnal->n;
  anal.BMR = pyAnal->BMR;
  anal.alpha = pyAnal->alpha;
  anal.prior_cols = pyAnal->prior_cols;
  pyRes->validResult = false;
  bool validModel = false;  // is a valid model selected for multitumor to run
  for (int dataset = 0; dataset < pyAnal->ndatasets; dataset++) {
    int selIndex = pyRes->selectedModelIndex[dataset];
    if (selIndex >= 0) {
      validModel = true;
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
      res.selectedModelIndex.push_back(0);  // selected index is always zero since only one model
      modRes = pyRes->models[dataset][selIndex];
      modResGroup.push_back(modRes);
      res.models.push_back(modResGroup);
    }
  }
  anal.ndatasets = anal.nmodels.size();
  res.ndatasets = anal.ndatasets;

  if (!validModel) {
    std::cout << "No valid models available for multitumor analysis" << std::endl;
    return;
  }

  // run MSCombo
  runMultitumorModel(&anal, &res);

  pyRes->BMD = res.BMD;
  pyRes->BMDL = res.BMDL;
  pyRes->BMDU = res.BMDU;
  pyRes->slopeFactor = res.slopeFactor;
  pyRes->combined_LL = res.combined_LL;
  pyRes->combined_LL_const = res.combined_LL_const;
  pyRes->validResult = true;
}

void selectMultitumorModel(
    struct python_multitumor_analysis *pyAnal, struct python_multitumor_result *pyRes
) {
  if (pyRes->selectedModelIndex.size() > 0) {
    pyRes->selectedModelIndex.clear();
  }

  for (int i = 0; i < pyAnal->ndatasets; i++) {
    if (pyAnal->degree[i] == 0) {
      int selectedIndex = selectBestMultitumorModel(pyAnal->models[i], pyRes->models[i]);
      pyRes->selectedModelIndex.push_back(selectedIndex);
    } else {
      // in this case the only the user selected model was run
      pyRes->selectedModelIndex.push_back(0);
    }
  }
}

int selectBestMultitumorModel(
    std::vector<python_dichotomous_analysis> &analModels,
    std::vector<python_dichotomous_model_result> &resModels
) {
  double bmdDoseSRLimit = 2.0;      // scaled residual at BMD dose < value
  double controlDoseSRLimit = 2.0;  // scaled residual at control dose < value
  double pValueMin = 0.05;          // p-value > value

  std::vector<bool> modelHitBoundary;  // checks whether specific model hit boundary
  std::vector<bool> adequateFit;
  bool anyHitBoundary = false;
  bool anyAdequateFit = false;

  for (int i = 0; i < analModels.size(); i++) {
    // check parameter boundaries
    modelHitBoundary.push_back(false);
    std::vector<bool> bounded = resModels[i].bmdsRes.bounded;
    for (int j = 0; j < bounded.size(); j++) {
      if (bounded[j]) {
        modelHitBoundary[i] = true;
      }
    }
    // check adequate fit
    // 1. scaled residual at BMD dose < bmdDoseSRLimit
    // 2. scaled residual at control dose < controlDoseSRLimit
    // 3. If 1 & 2 are met, p-value > pValueMin
    // Note: these values are all currently hardcoded.  Need to account for/allow user changes
    adequateFit.push_back(false);
    if (resModels[i].getSRAtDose(resModels[i].bmdsRes.BMD, analModels[i].doses) < bmdDoseSRLimit &&
        resModels[i].gof.residual[0] < controlDoseSRLimit && resModels[i].gof.p_value > pValueMin) {
      adequateFit[i] = true;
    }
  }

  // first only examine models up to k-2
  anyHitBoundary = false;
  anyAdequateFit = false;
  for (int i = 0; i < analModels.size() - 1; i++) {
    if (modelHitBoundary[i]) {
      anyHitBoundary = true;
    }
    if (adequateFit[i]) {
      anyAdequateFit = true;
    }
  }

  if (anyHitBoundary) {
    // at least one model had a parameter that hit a boundary
    // step 2 under SOP instructions
    if (adequateFit[0] && adequateFit[1]) {
      // both linear and quadratic models fit
      if (!modelHitBoundary[0] && !modelHitBoundary[1]) {
        // choose model with lowest AIC
        if (round_to(resModels[1].bmdsRes.AIC, BMDS_EPS) <
            round_to(resModels[0].bmdsRes.AIC, BMDS_EPS)) {
          return 1;
        } else {
          return 0;
        }
      } else {
        // choose model with lowest BMDL
        if (round_to(resModels[1].bmdsRes.BMDL, BMDS_EPS) <
            round_to(resModels[0].bmdsRes.BMDL, BMDS_EPS)) {
          return 1;
        } else {
          return 0;
        }
      }
    } else if (adequateFit[0]) {
      // only linear fits
      return 0;
    } else if (adequateFit[1]) {
      // only quadratic fits
      return 1;
    } else {
      // neither fit - refer to SWG
      return BMDS_MISSING;
    }
    // return 1;

  } else {
    // no models had a paramter hit the boundary
    // step 1 under SOP instructions
    if (anyAdequateFit) {
      // use lowest AIC model that adequately fits up to k-2 model (exclude k-1 model)
      std::vector<double> modelAIC;
      for (int i = 0; i < resModels.size() - 1; i++) {
        modelAIC.push_back(resModels[i].bmdsRes.AIC);
      }
      // now find the lowest AIC
      int index = std::distance(
          std::begin(modelAIC), std::min_element(std::begin(modelAIC), std::end(modelAIC))
      );
      return index;
    } else {
      // no adequateFit models.  Examine k-1 model
      if (adequateFit[analModels.size() - 1]) {
        return analModels.size() - 1;  // use k-1 model
      } else {
        return BMDS_MISSING;  // no model selected
      }
    }
  }
}

double calcSlopeFactor(double bmr, double bmdl) {
  if (bmdl > 0.0) {
    return bmr / bmdl;
  } else {
    return BMDS_MISSING;
  }
}

void BMDS_ENTRY_API __stdcall runMultitumorModel(
    struct python_multitumor_analysis *pyAnal, struct python_multitumor_result *pyRes
) {
  // calculate slopeFactor for all individual models
  for (int i = 0; i < pyAnal->ndatasets; i++) {
    for (int j = 0; j < pyAnal->nmodels[i]; j++) {
      double bmr = pyAnal->models[i][j].BMR;
      pyRes->models[i][j].bmdsRes.setSlopeFactor(bmr);
    }
  }

  // calculate combined LL & constant
  pyRes->combined_LL = 0.0;
  pyRes->combined_LL_const = 0.0;
  for (int i = 0; i < pyAnal->ndatasets; i++) {
    int selIndex = pyRes->selectedModelIndex[i];
    pyRes->combined_LL += pyRes->models[i][selIndex].max;
    pyRes->combined_LL_const +=
        LogLik_Constant(pyAnal->models[i][selIndex].Y, pyAnal->models[i][selIndex].n_group);
  }
  pyRes->combined_LL = pyRes->combined_LL * -1;
  pyRes->combined_LL_const = pyRes->combined_LL_const * -1;

  Multistage_ComboBMD(pyAnal, pyRes);

  pyRes->setSlopeFactor(pyAnal->BMR);
}

double PROBABILITY_INRANGE(double ex) {
  if (ex <= 0.00000001) {
    ex = 1.0e-7;
  }
  if (ex >= 0.9999999) {
    ex = 0.9999999;
  }
  return ex;
}

void BMDS_ENTRY_API __stdcall pythonBMDSNested(
    struct python_nested_analysis *pyAnal, struct python_nested_result *pyRes
) {
  //////////////////////////////
  // Parm vector indexes
  // 0 = alpha (intercept or background)
  // 1 = beta
  // 2 = THETA1
  // 3 = THETA2
  // 4 = rho
  // 5 = PHISTART - first PHI parameter
  //
  int BGINDEX = 0;       // location of background parameter in parms vector
  int PHISTART = 5;      // location of first PHI parameter in parms vector
  int THETA1_INDEX = 2;  // location of THETA1 parameter in parms vector
  int THETA2_INDEX = 3;  // location of THETA2 parameter in parms vector
  pyRes->validResult = false;
  pyRes->bmd = BMDS_MISSING;
  pyRes->bmdsRes.BMD = BMDS_MISSING;
  pyRes->bmdsRes.BMDL = BMDS_MISSING;
  pyRes->bmdsRes.BMDU = BMDS_MISSING;

  // set seed from time clock if default seed=0 is specified
  if (pyAnal->seed == BMDS_MISSING) {
    pyAnal->seed = time(NULL);
  }

  int Nobs = pyAnal->doses.size();
  std::vector<double> Xi(Nobs);  // doses (independent data var)
  std::vector<double> Yp(Nobs);  // incidence (positive dependent var)
  std::vector<double> Yn(Nobs);  // negative dependent var
  std::vector<double> Ls(Nobs);  // litter size
  std::vector<double> SR(Nobs);  // scaled residual
  std::vector<double> Lsc(Nobs);
  std::vector<int> Xg(Nobs, BMDS_MISSING);  // markings for group number

  Xi = pyAnal->doses;
  Yp = pyAnal->incidence;
  for (int i = 0; i < Nobs; i++) {
    Yn[i] = pyAnal->litterSize[i] - pyAnal->incidence[i];
  }
  Ls = pyAnal->litterSize;
  Lsc = pyAnal->lsc;

  double junk1 = Xi[0];
  int junk = 1;
  for (int i = 0; i < Nobs; i++) {
    if (Xi[i] == junk1) {
      Xg[i] = junk;
    } else {
      Xg[i] = ++junk;
      junk1 = Xi[i];
    }
  }
  int ngrp = junk;
  // compute the group size
  std::vector<int> grpSize(ngrp);
  for (int i = 0; i < Nobs; i++) {
    grpSize[Xg[i] - 1] += 1;
  }
  SortNestedData(grpSize, Xi, Ls, Yp, Yn, Lsc, false);  // Sort by Xi->Ls->Yp

  // push last group size;
  pyRes->nparms = 5 + ngrp;
  std::vector<bool> Spec(pyRes->nparms, false);

  // create and initialize vectors needed
  pyRes->parms.resize(pyRes->nparms);
  for (int i = 0; i < pyRes->nparms; i++) {
    pyRes->parms[i] = BMDS_MISSING;
  }
  pyRes->bmdsRes.stdErr.resize(pyRes->nparms, BMDS_MISSING);
  pyRes->bmdsRes.bounded.resize(pyRes->nparms, false);
  pyRes->litter.dose.resize(Nobs, BMDS_MISSING);
  pyRes->litter.LSC.resize(Nobs, BMDS_MISSING);
  pyRes->litter.estProb.resize(Nobs, BMDS_MISSING);
  pyRes->litter.litterSize.resize(Nobs, BMDS_MISSING);
  pyRes->litter.expected.resize(Nobs, BMDS_MISSING);
  pyRes->litter.observed.resize(Nobs, static_cast<int>(BMDS_MISSING));
  pyRes->litter.SR.resize(Nobs, BMDS_MISSING);
  pyRes->boot.pVal.resize(pyAnal->numBootRuns + 1, BMDS_MISSING);
  pyRes->boot.perc50.resize(pyAnal->numBootRuns + 1, BMDS_MISSING);
  pyRes->boot.perc90.resize(pyAnal->numBootRuns + 1, BMDS_MISSING);
  pyRes->boot.perc95.resize(pyAnal->numBootRuns + 1, BMDS_MISSING);
  pyRes->boot.perc99.resize(pyAnal->numBootRuns + 1, BMDS_MISSING);
  pyRes->reduced.dose.resize(ngrp, BMDS_MISSING);
  pyRes->reduced.propAffect.resize(ngrp, BMDS_MISSING);
  pyRes->reduced.lowerConf.resize(ngrp, BMDS_MISSING);
  pyRes->reduced.upperConf.resize(ngrp, BMDS_MISSING);

  int knownParms = 0;
  // handle specified parms
  if (!pyAnal->estBackground) {
    // set background to zero
    pyRes->parms[BGINDEX] = 0;
    Spec[BGINDEX] = true;
    knownParms++;
  }
  if (pyAnal->ILC_type == 0) {
    // set phi parm values to zero
    for (int i = PHISTART; i < pyRes->nparms; i++) {
      pyRes->parms[i] = 0.0;
      Spec[i] = true;
    }
    knownParms += ngrp;
  }
  if (pyAnal->LSC_type == 0) {
    // set theta parm values to zero
    for (int i = THETA1_INDEX; i <= THETA2_INDEX; i++) {
      pyRes->parms[i] = 0.0;
      Spec[i] = true;
    }
    knownParms += 2;  // theta1 and theta2
  }

  if (Nobs < pyRes->nparms - knownParms) {
    std::cout << "Error: Fewer observations " << Nobs << " than estimated parameters "
              << pyRes->nparms - knownParms << std::endl;
    return;
  }

  // make sure all litter sizes > 0
  for (int i = 0; i < Nobs; i++) {
    if (pyAnal->litterSize[i] <= 0.0) {
      std::cout << "Error: All litter sizes must be greater than zero." << std::endl;
      return;
    }
  }

  double EPS = 3.0e-8;
  int iter = 0;
  double xlk = 0.0;
  ////////////////////////////////////////////////////
  // START NLOGIST_FIT
  // //////////////////////////////////////////////////

  double sum1 = 0.0;
  double nij = 0.0;
  int index = 0;
  while (Xg[index] == 1) {
    sum1 += Lsc[index];
    nij += 1.0;
    index++;
  }
  double smean1 = sum1 / nij;  // the average lsc in group 1

  sum1 = 0.0;
  nij = 0.0;
  double smax = Lsc[0];
  double smin = Lsc[0];
  double xmax = 0.0;
  double x = 0.0;
  for (int i = 0; i < Nobs; i++) {
    x = Xi[i];
    sum1 += Lsc[i];
    nij += 1.0;
    if (Ls[i] > smax) smax = Lsc[i];
    if (Ls[i] < smin) smin = Lsc[i];
    if (x > xmax) xmax = x;
  }
  double smean = sum1 / nij;  // overall average lsc

  double sijfixed =
      smean;  // default to overall mean.  Not used if pyAnal->LSCType = 0 (do not use LSC)
  if (pyAnal->LSC_type == 2) {  // control group mean
    sijfixed = smean1;
  }
  pyRes->fixedLSC = sijfixed;

  // compute default starting values for all unfixed parameters
  // revisit if we need to fix any parameters

  // Get the starting values in stages
  // First, get starting values for alpha, beta, and rho

  std::vector<double> tmpXi(Nobs);  // doses (independent data var)
  std::vector<double> tmpYp(Nobs);  // incidence (positive dependent var)
  std::vector<double> tmpYn(Nobs);  // negative dependent var
  std::vector<double> tmpy(2);
  Eigen::Matrix2d tmpvcv(2, 2);

  index = 0;
  // convert nested data to dichotomous data
  for (int j = 0; j < ngrp; j++) {
    tmpYn[j] = 0;
    tmpYp[j] = 0;
    while (index < Nobs && Xg[index] == (j + 1)
    ) {  // note - index is 1 behind group number (index=0 for group 1)
      tmpYn[j] += Yn[index];
      tmpYp[j] += Yp[index];
      tmpXi[j] = Xi[index];
      index++;
    }
  }

  // compute initial estimates
  double ymin = 1.0;
  double W = 0.0;
  for (int i = 0; i < ngrp; i++) {
    W = tmpYp[i] / (tmpYp[i] + tmpYn[i]);
    probability_inrange(&W);
    if (W < ymin) ymin = W;
  }

  for (int i = 0; i < ngrp; i++) {
    x = tmpXi[i];
    W = log((tmpYp[i] - ymin * (tmpYp[i] + tmpYn[i]) + 0.5) / (tmpYn[i] + 0.5));
    if (x <= CloseToZero) {
      x = Log_zero;
    } else {
      x = log(x);
    }
    tmpvcv(0, 0) += 1.0;
    tmpvcv(0, 1) += x;
    tmpvcv(1, 1) += x * x;
    tmpy[0] += W;
    tmpy[1] += W * x;
  }
  tmpvcv(1, 0) = tmpvcv(0, 1);
  pyRes->parms[0] = ymin + 0.001;  // in case ymin=0

  if (!Spec[0] && !Spec[1]) {
    if (tmpvcv.determinant() > 0) {
      Eigen::Matrix2d invtmpvcv = tmpvcv.inverse();
      pyRes->parms[1] = invtmpvcv(0, 0) * tmpy[0] + invtmpvcv(0, 1) * tmpy[1];
      pyRes->parms[4] = invtmpvcv(1, 0) * tmpy[0] + invtmpvcv(1, 1) * tmpy[1];
    } else {
      pyRes->parms[1] = -1.0;
      pyRes->parms[4] = 1.001;
    }
  }

  if (pyRes->parms[4] < pyAnal->prior[4]) {
    pyRes->parms[4] = 1.0001;
  }

  // setup log-logistic model
  for (int i = 5; i < pyRes->nparms; i++) {
    pyRes->parms[i] = 0.0;
  }
  pyRes->parms[2] = pyRes->parms[3] = 0.0;

  // simple non-nested log-logistic model

  // This first pass has theta1 p[2] and theta2 p[3] fixed to zero
  struct nestedObjData objData;
  objData.Ls = Ls;
  objData.Xi = Xi;
  objData.Xg = Xg;
  objData.Yp = Yp;
  objData.Yn = Yn;
  objData.Lsc = Lsc;
  objData.ngrp = ngrp;
  objData.LSC_type = pyAnal->LSC_type;
  objData.ILC_type = pyAnal->ILC_type;
  objData.smax = smax;
  objData.smin = smin;
  objData.isBMDL = false;
  objData.sijfixed = sijfixed;
  objData.riskType = pyAnal->BMD_type;
  objData.BMR = pyAnal->BMR;
  objData.tol = 1e-8;
  objData.prior = pyAnal->prior;
  objData.Spec = Spec;
  objData.Spec[THETA1_INDEX] = true;
  objData.Spec[THETA2_INDEX] = true;
  for (int i = PHISTART; i < pyRes->nparms; i++) {
    objData.Spec[i] = true;
  }

  std::vector<double> pBak = pyRes->parms;

  double retVal = opt_nlogistic(pyRes->parms, &objData);

  if (retVal < 0) return;  // return with validResult = false

  // Now alpha p[0], beta p[1], and rho[4] contain starting estimates
  // Second, get initial values for Phi's
  int count = 0;
  for (int i = PHISTART; i < pyRes->nparms; i++) {
    if (Spec[i]) count++;
  }
  // Currently do not allow specification of Phi values, so this is always true
  if (count < ngrp && pyAnal->ILC_type != 0) {
    // leave theta1 and theta2 = 0.0;
    for (int i = 5; i < pyRes->nparms; i++) {
      // set back to original parms
      // if (SpBak[i] == 0){
      // pyRes->parms[i] = 0.01;
      pyRes->parms[i] = pBak[i];
      //}
    }
    for (int i = 5; i < pyRes->nparms; i++) {
      pyRes->parms[i] = pyRes->parms[i] / (1 - pyRes->parms[i]);  // Phi --> Psi
    }

    objData.Spec = Spec;
    objData.Spec[THETA1_INDEX] = true;
    objData.Spec[THETA2_INDEX] = true;
    // pass 2 has thetas set to zero

    // allow to continue even if retVal < 0
    retVal = opt_nlogistic(pyRes->parms, &objData);

    // Transform parameters to "external" form
    for (int i = 5; i < pyRes->nparms; i++) {
      pyRes->parms[i] = pyRes->parms[i] / (1 + pyRes->parms[i]);  // Psi --> Phi
    }

    // When theta1 ==0, internal and external forms for alpha are the same,
    // so we do not need to back transform p[0]
  }

  // Finally, get starting values for Thetas
  for (int i = 2; i <= 3; i++) {
    // if not set
    pyRes->parms[i] = 0.0;
  }

  // Fit the model

  // Transform the parameters to the "internal" form
  for (int i = 5; i < pyRes->nparms; i++) {
    pyRes->parms[i] = pyRes->parms[i] / (1 - pyRes->parms[i]);  // Phi --> Psi
  }

  double junk3;
  if (Spec[0] == Spec[2]) {
    junk1 = pyRes->parms[0];
    junk3 = pyRes->parms[2];
    pyRes->parms[0] = junk1 + smin * junk3;
    pyRes->parms[2] = junk1 + smax * junk3;
  }

  // ML fit and return log-likelihood value
  objData.Spec = Spec;

  for (int i = 0; i < pyRes->nparms; i++) {
    if (Spec[i]) {
      pyRes->parms[i] = 0.0;
    }
  }

  if (pyRes->parms[4] <= pyAnal->prior[4]) {
    pyRes->parms[4] = 1.0001;
  }
  objData.tol = 1e-12;
  retVal = opt_nlogistic(pyRes->parms, &objData);
  objData.BMD_lk = objData.xlk;  // save BMD likelihood
  pyRes->LL = objData.xlk;

  // Transform the parameters to the "external" form
  for (int i = 5; i < pyRes->nparms; i++) {
    pyRes->parms[i] = pyRes->parms[i] / (1 + pyRes->parms[i]);
  }

  double sdif = smax - smin;
  if (Spec[0] == Spec[2]) {
    junk1 = pyRes->parms[0];
    junk3 = pyRes->parms[2];
    pyRes->parms[0] = (smax * junk1 - smin * junk3) / sdif;
    pyRes->parms[2] = (junk3 - junk1) / sdif;
  }

  // TODO:  check retVal to make sure convergence was achieved
  if (retVal > 0) {
    pyRes->validResult = true;
    pyRes->bmdsRes.validResult = true;
  }

  // compute Hessian
  Eigen::MatrixXd vcv(pyRes->nparms, pyRes->nparms);

  std::vector<std::vector<double>> vcv_tmp(pyRes->nparms, std::vector<double>(pyRes->nparms));
  Nlogist_vcv(pyRes->parms, pyRes->bmdsRes.bounded, &objData, vcv_tmp);

  for (int i = 0; i < pyRes->nparms; i++) {
    vcv.row(i) = Eigen::VectorXd::Map(&vcv_tmp[i][0], vcv_tmp[i].size());
  }

  std::vector<double> lowerBound(pyAnal->prior.begin(), pyAnal->prior.begin() + pyRes->nparms);
  std::vector<double> upperBound(pyAnal->prior.begin() + pyRes->nparms, pyAnal->prior.end());

  // TODO:  make sure this catches all parms in Spec vector
  int numBounded = checkForBoundedParms(
      pyRes->nparms, &pyRes->parms[0], &lowerBound[0], &upperBound[0], &pyRes->bmdsRes
  );

  // remove rows and columns of vcv that correspond to bounded parms
  for (int i = 0; i < pyRes->nparms; i++) {
    if (pyRes->bmdsRes.bounded[i]) {
      vcv.row(i).setZero();
      vcv.col(i).setZero();
    }
  }

  Eigen::MatrixXd red_vcv(pyRes->nparms - numBounded, pyRes->nparms - numBounded);
  int jvar = 0;
  int ivar = 0;
  for (int i = 0; i < vcv.cols(); i++) {
    if (!pyRes->bmdsRes.bounded[i]) {
      for (int j = 0; j < vcv.rows(); j++) {
        if (!pyRes->bmdsRes.bounded[j]) {
          red_vcv(ivar, jvar) = vcv(i, j);
          jvar++;
        }
      }
      ivar++;
      jvar = 0;
    }
  }

  Eigen::MatrixXd inv_vcv(pyRes->nparms, pyRes->nparms);
  if (red_vcv.determinant() > 0) {
    inv_vcv = red_vcv.inverse();
    int j = 0;
    for (int i = 0; i < pyRes->nparms; i++) {
      if (!pyRes->bmdsRes.bounded[i] && inv_vcv(j, j) > 0) {
        pyRes->bmdsRes.stdErr[i] = sqrt(inv_vcv(j, j));
        j++;
      } else {
        if (j < (inv_vcv.rows() - 1) && inv_vcv(j, j) < 0) {
          j++;
        }
        pyRes->bmdsRes.stdErr[i] = BMDS_MISSING;
      }
    }
  }

  // compute and output the analysis of deviance

  struct nested_AOD fullAnova;
  struct nested_AOD fittedAnova;
  struct nested_AOD redAnova;

  // compute init_lkf for full model and init_lkr for reduced model
  double lkf = 0.0;   // full model likelihood
  double pdep = 0.0;  // positive response
  double ndep = 0.0;  // negative response

  for (int i = 0; i < Nobs; i++) {
    pdep += Yp[i];
    ndep += Yn[i];
    W = Yp[i] / (Yp[i] + Yn[i]);
    if (W > 0) {
      lkf += Yp[i] * log(W);
    }
    if (W < 1) {
      lkf += Yn[i] * log(1 - W);
    }
  }
  W = pdep / (pdep + ndep);
  double lkr = pdep * log(W) + ndep * log(1 - W);  // reduced model likelihood

  double dev = 2 * (lkf - pyRes->LL);
  int df = pyRes->nparms - numBounded;
  double pv;
  if (dev < 0.0) {
    pv = BMDS_MISSING;
  } else {
    if (df > 0) {
      pv = CHISQ(dev, df);
    } else {
      pv = BMDS_MISSING;
    }
  }

  int numSpec = 0;
  for (int i = 0; i < Spec.size(); i++) {
    if (Spec[i]) numSpec++;
  }

  fullAnova.df = pyRes->nparms - numSpec;
  fullAnova.LL = lkf;

  fittedAnova.LL = pyRes->LL;
  fittedAnova.dev = dev;
  fittedAnova.df = Nobs - df;
  fittedAnova.pv = pv;

  if (Nobs > 1) {
    dev = 2 * (lkf - lkr);
    if (dev <= 0) {
      pv = BMDS_MISSING;
    } else {
      pv = CHISQ(dev, Nobs - 1);
    }
  } else {
    pv = BMDS_MISSING;
  }

  redAnova.LL = lkr;
  redAnova.dev = dev;
  redAnova.df = Nobs - 1;
  redAnova.pv = pv;

  // Vector containing
  std::vector<nested_AOD> anovaVec;
  anovaVec.push_back(fullAnova);
  anovaVec.push_back(fittedAnova);
  anovaVec.push_back(redAnova);

  pyRes->bmdsRes.AIC = -2 * fittedAnova.LL + 2 * (1.0 + redAnova.df - fittedAnova.df);

  pyRes->model_df = fittedAnova.df;

  // Generate litter data results
  Nlogist_GOF(pyRes->parms, &objData, &pyRes->litter, grpSize);
  pyRes->bmdsRes.chisq = pyRes->litter.chiSq;

  std::vector<double> GXi(ngrp);
  for (int j = 0; j < Nobs; j++) {
    for (int i = 0; i < ngrp; i++) {
      if (Xg[j] == i + 1) {
        GXi[i] = Xi[j];
      }
    }
  }

  objData.GXi = GXi;

  Nlogist_reduced(pyAnal->alpha, &objData, &pyRes->reduced);

  // Compute BMD
  Nlogist_BMD(pyAnal, pyRes, smin, smax, sijfixed, xmax, &objData);

  Nlogist_SRoI(&objData, &pyRes->srData, pyRes->litter.SR, grpSize, pyRes->bmd);

  Nlogist_Bootstrap(&objData, pyRes, pyAnal->seed, pyAnal->iterations, pyAnal->numBootRuns);
}

void Nlogist_Bootstrap(
    struct nestedObjData *objData, struct python_nested_result *pyRes, long seed, int iterations,
    int BSLoops
) {
  int ngrp = objData->ngrp;
  std::vector<double> Yp = objData->Yp;
  std::vector<double> Yn = objData->Yn;
  std::vector<int> Xg = objData->Xg;
  int Nobs = Yp.size();
  std::vector<double> grpSize(ngrp, 0);
  std::vector<double> phi(ngrp);
  std::vector<double> Ep(Nobs);
  std::vector<double> Ysum(Nobs);
  std::vector<double> Ypp(Nobs);
  std::vector<double> var(Nobs);

  // bootstrap variables
  gsl_rng *r;                                                  // random value for beta distribution
  const gsl_rng_type *type;                                    // type for beta distribution
  double urand;                                                /* uniform random value [0,1]  */
  int SRsqCount;                                               // tracks SR^2 >= SR^2 original
  std::vector<double> percentiles = {0.50, 0.90, 0.95, 0.99};  // percentile values to calculate
  std::vector<std::vector<double>> SR_newsqsum(BSLoops, std::vector<double>(iterations, 0));
  std::vector<double> Yp_new(Nobs);  // pseudo-observed values
  std::vector<double> SR_new(Nobs);  // SR values based on new variables
  double cutpoint;                   // cut-off probability for bootstrapping
  double a, b;                       // values for beta distribution function
  double x, cum_beta, cum_beta_comp, cum_beta_comp_low,
      cum_beta_comp_hi;  // cumulative beta variables for NaN method
  std::vector<double> pv(BSLoops);
  double cum_beta_step = 0.00001;  // step size for cumulative beta table
  double eps = 10e-8;              // accuracy to test for est.prob = 0;
  std::vector<double> pctTemp(BSLoops + 1);
  int perloc;

  // compute the # of obs in each group
  for (int i = 0; i < Nobs; i++) grpSize[Xg[i] - 1] += 1;
  // assumes PHI_START = 5
  for (int i = 0; i < ngrp; i++) {
    phi[i] = pyRes->parms[5 + i];
  }

  // compute the estimated probability and expected obs
  Nlogist_Predict(pyRes->parms, objData, Ep);

  for (int i = 0; i < Nobs; i++) {
    Ysum[i] = Yp[i] + Yn[i];
    Ypp[i] = Ep[i] * Ysum[i];
    var[i] = Ysum[i] * Ep[i] * (1 - Ep[i]) * (1.0 + (Ysum[i] - 1.0) * phi[Xg[i] - 1]);
  }

  double gfit = pyRes->litter.chiSq;
  int df = pyRes->model_df;

  // bootstrap method for nested models

  // set up GSL random number generator
  gsl_rng_env_setup();
  type = gsl_rng_default;
  r = gsl_rng_alloc(type);
  gsl_rng_set(r, seed);

  for (int l = 0; l < BSLoops; l++) {
    SRsqCount = 0;
    for (int i = 0; i < iterations; i++) {
      for (int j = 0; j < Nobs; j++) {  // loop over each line in litter data
        Yp_new[j] = 0;                  // reset pseudo-observed value to zero

        // find cut-off probability (cutpoint)
        if ((phi[Xg[j] - 1] == 0) || (phi[Xg[j] - 1] == 1)) {
          cutpoint = Ep[j];
        } else {
          // find cut-off probability from beta distribution
          // be sure not to evaluate a,b and beta distribution if Ep[j] = 0.  This will skip next if
          // statement also: if (isnan(cutpoint))
          if (Ep[j] <= eps) {
            cutpoint = 0.0;
          } else {
            double tempvalue =
                Ep[j] * (1.0 - Ep[j]) / (Ep[j] * (1.0 - Ep[j]) * phi[Xg[j] - 1]) - 1.0;
            a = Ep[j] * tempvalue;
            b = (1.0 - Ep[j]) * tempvalue;
            cutpoint = gsl_ran_beta(r, a, b);
          }
        }
        if (isnan(cutpoint)) {
          cum_beta_comp_low = gsl_cdf_beta_P(cum_beta_step, a, b);  // low end of cumulative scale
          cum_beta_comp_hi =
              gsl_cdf_beta_P(1.0 - cum_beta_step, a, b);  // high end of cumulative scale
          cum_beta = gsl_rng_uniform(r);  // random value to compare to cumulative scale

          if (cum_beta < cum_beta_comp_low) {  // if random value is less than low end
            cutpoint = 0;
          } else if (cum_beta >= cum_beta_comp_hi) {  // if random value is greater than high end
            cutpoint = 1.0 - cum_beta_step;
          } else if (fabs(cum_beta_comp_low - cum_beta) <
                     fabs(cum_beta_comp_hi - cum_beta)) {  // if random value is closer to low end
            cum_beta_comp = 0;
            x = 0;
            while (x <= 1.0 && cum_beta > cum_beta_comp) {
              cutpoint = x;
              x += cum_beta_step;
              cum_beta_comp = gsl_cdf_beta_P(x, a, b);
            }
          } else {  // random value is closer to high end
            cum_beta_comp = 1.0;
            x = 1;
            cutpoint = 1.0;
            while (x >= 0.0 && cum_beta < cum_beta_comp) {
              x -= cum_beta_step;
              cutpoint = x;
              cum_beta_comp = gsl_cdf_beta_P(x, a, b);
            }
          }
        }
        // calculate pseudo-observed value
        if (Ep[j] <= eps) {  // if est prob = 0, then bootstrap must produce zero responses for this
                             // observation
          Yp_new[j] = 0;
        } else if (phi[Xg[j] - 1] != 1) {
          for (int k = 0; k < Ysum[j]; k++) {
            urand = gsl_rng_uniform(r);  // generate uniform random number
            if (urand <= cutpoint) Yp_new[j]++;
          }
        } else {  // phi = 1 method
          urand = gsl_rng_uniform(r);
          if (urand <= cutpoint) {
            Yp_new[j] = Ysum[j];
          } else {
            Yp_new[j] = 0;
          }
        }

        // calculate scaled residual based on pseudo-observation Yp_new
        SR_new[j] = var[j] > 0 ? (Yp_new[j] - Ypp[j]) / sqrt(var[j]) : 0;
        SR_newsqsum[l][i] += SR_new[j] * SR_new[j];

      }  // end Nobs loop

      if (SR_newsqsum[l][i] >= gfit) SRsqCount++;

    }  // end iterations loop

    pv[l] = (double)SRsqCount / iterations;  // compute p-value for each BS loop

  }  // end BSLoops loop

  // calculate chi-square percentiles for each BS loop
  double pavg = 0;

  for (int l = 0; l < BSLoops; l++) {
    pyRes->boot.pVal[l] = pv[l];
    // compute percentile values of sum of SR^2
    sort(SR_newsqsum[l].begin(), SR_newsqsum[l].end());
    for (int k = 0; k < percentiles.size(); k++) {
      double temp = percentiles[k] * iterations;
      if ((ceilf(temp) == temp) && (floorf(temp) == temp)) {
        perloc = (int)temp;
        pctTemp[k] = SR_newsqsum[l][perloc];
      } else {
        perloc = (int)ceilf(temp);
        pctTemp[k] = (SR_newsqsum[l][perloc] + SR_newsqsum[l][perloc + 1]) / 2.0;
      }
    }
    pyRes->boot.perc50[l] = pctTemp[0];
    pyRes->boot.perc90[l] = pctTemp[1];
    pyRes->boot.perc95[l] = pctTemp[2];
    pyRes->boot.perc99[l] = pctTemp[3];
    pavg += pv[l];
  }

  // combined p-value
  pavg /= BSLoops;
  pyRes->boot.pVal[BSLoops] = pavg;
  pyRes->combPVal = pavg;

  // combined chi-square percentiles
  // 1st flatten the 2d vector
  std::vector<double> combSR(BSLoops * iterations);
  int count = 0;
  for (int i = 0; i < BSLoops; i++) {
    for (int j = 0; j < iterations; j++) {
      combSR[count] = SR_newsqsum[i][j];
      count++;
    }
  }
  sort(combSR.begin(), combSR.end());
  for (int k = 0; k < percentiles.size(); k++) {
    double temp = percentiles[k] * BSLoops * iterations;
    if ((ceilf(temp) == temp) && (floorf(temp) == temp)) {
      perloc = (int)temp;
      pctTemp[k] = combSR[perloc - 1];  // CHECK THIS!!!
    } else {
      perloc = (int)ceilf(temp);
      pctTemp[k] = (combSR[perloc - 1] + combSR[perloc]) / 2.0;
    }
  }
  pyRes->boot.perc50[BSLoops] = pctTemp[0];
  pyRes->boot.perc90[BSLoops] = pctTemp[1];
  pyRes->boot.perc95[BSLoops] = pctTemp[2];
  pyRes->boot.perc99[BSLoops] = pctTemp[3];

  gsl_rng_free(r);  // free memory from random generator
}

void validatePositiveInput(double &val) {
  if (val < 0) {
    val = BMDS_EPS;
  }
}

void Nlogist_SRoI(
    struct nestedObjData *objData, struct nestedSRData *srData, std::vector<double> &SR,
    const std::vector<int> &grpSize, double bmd
) {
  int locDose = 0;
  int locLSC = 0;
  double closeDose = 0.0;
  double closeLSC = 0.0;
  int litSR = 1;
  double idiff;

  int ngrp = objData->ngrp;
  std::vector<double> Xi = objData->Xi;
  std::vector<int> Xg = objData->Xg;
  std::vector<double> Lsc = objData->Lsc;
  std::vector<double> Yp = objData->Yp;
  std::vector<double> Yn = objData->Yn;
  std::vector<double> Ls = objData->Ls;
  std::vector<double> GXi = objData->GXi;

  // Need to sort data by Lsc
  SortNestedData(grpSize, Xi, Ls, Yp, Yn, Lsc, true);

  double meanLSC = objData->sijfixed;
  int Nobs = Xi.size();

  // choose dose group closest to BMD
  double diff = DBL_MAX;
  for (int i = 0; i < ngrp; i++) {
    idiff = fabs(bmd - GXi[i]);
    if (idiff < diff) {
      diff = idiff;
      locDose = i;
      closeDose = GXi[i];
    }
  }

  // choose LSC closest to mean LSC
  diff = DBL_MAX;
  for (int i = 0; i < Nobs; i++) {
    if (Xi[i] == closeDose) {
      idiff = fabs(Lsc[i] - meanLSC);
      if (idiff == diff) litSR++;
      if (idiff < diff) {
        litSR = 1;
        diff = idiff;
        closeLSC = Lsc[i];
        locLSC = i;
      }
    }
  }

  // calculate max, min, average SRoI and |SRoI|
  srData->maxSR = SR[locLSC];
  srData->minSR = SR[locLSC];
  srData->avgSR = SR[locLSC];
  srData->maxAbsSR = fabs(SR[locLSC]);
  srData->minAbsSR = fabs(SR[locLSC]);
  srData->avgAbsSR = fabs(SR[locLSC]);

  if (litSR != 1) {
    srData->avgSR = 0;
    srData->avgAbsSR = 0;
    for (int i = locLSC; i <= locLSC + litSR - 1;
         i++) {  // relies on dose groups being sorted by LSC
      if (SR[i] > srData->maxSR) srData->maxSR = SR[i];
      if (SR[i] < srData->minSR) srData->minSR = SR[i];
      if (fabs(SR[i]) > srData->maxAbsSR) srData->maxAbsSR = fabs(SR[i]);
      if (fabs(SR[i]) < srData->minAbsSR) srData->minAbsSR = fabs(SR[i]);
      srData->avgSR += SR[i];
      srData->avgAbsSR += fabs(SR[i]);
    }
    srData->avgSR /= litSR;
    srData->avgAbsSR /= litSR;
  }
}

void Nlogist_reduced(
    double alpha, struct nestedObjData *objData, struct nestedReducedData *redData
) {
  int ngrp = objData->ngrp;
  std::vector<double> num(ngrp);
  std::vector<double> den(ngrp);
  raoscott(objData, num, den);

  redData->dose = objData->GXi;
  // Quantal_CI requires number pos and number neg
  // after this, den contains the number unaffected (transformed)
  for (int i = 0; i < ngrp; i++) den[i] -= num[i];
  double conf = 1 - alpha;
  Quantal_CI(num, den, conf, redData);
}

void Quantal_CI(
    std::vector<double> &Yp, std::vector<double> &Yn, double conf, struct nestedReducedData *redData
) {
  double n, phat;
  double p = 1.0 - (1.0 - conf) / 2.0;
  double sigma = 1.0;
  double cc = gsl_cdf_gaussian_Pinv(p, sigma);
  double cc2 = cc * cc;
  int Nobs = Yp.size();

  for (int i = 0; i < Nobs; i++) {
    n = Yp[i] + Yn[i];
    redData->propAffect[i] = phat = Yp[i] / n;
    redData->lowerConf[i] =
        (2 * Yp[i] + cc2 - 1.0) - cc * sqrt(cc2 - (2.0 + 1.0 / n) + 4.0 * phat * (Yn[i] + 1.0));
    redData->lowerConf[i] /= 2.0 * (n + cc2);
    redData->upperConf[i] =
        (2 * Yp[i] + cc2 + 1.0) + cc * sqrt(cc2 + (2.0 - 1.0 / n) + 4.0 * phat * (Yn[i] - 1.0));
    redData->upperConf[i] /= 2.0 * (n + cc2);
  }
}

void raoscott(struct nestedObjData *objData, std::vector<double> &num, std::vector<double> &den) {
  int ngrp = objData->ngrp;
  int Nobs = objData->Xi.size();
  std::vector<int> Xg = objData->Xg;
  std::vector<double> Yp = objData->Yp;
  std::vector<double> Yn = objData->Yn;
  double sumnum, sumden, sumsq, phat, rij, vi, di;
  int grpsize;

  for (int j = 1; j <= ngrp; j++) {
    sumnum = sumden = 0.0;
    for (int i = 0; i < Nobs; i++) {
      if (j == Xg[i]) {
        sumnum += Yp[i];
        sumden += (Yp[i] + Yn[i]);
      }
    }
    phat = sumnum / sumden;
    sumsq = 0;
    grpsize = 0;
    for (int i = 0; i < Nobs; i++) {
      if (j == Xg[i]) {
        grpsize++;
        rij = Yp[i] - phat * (Yp[i] + Yn[i]);
        sumsq += rij * rij;
      }
    }
    vi = grpsize * sumsq / ((grpsize - 1) * sumden * sumden);
    if (phat > 0 && phat < 1) {
      di = sumden * vi / (phat * (1.0 - phat));
    } else {
      di = 1.0;
    }
    num[j - 1] = sumnum / di;
    den[j - 1] = sumden / di;
  }
}

void Nlogist_GOF(
    const std::vector<double> &parms, struct nestedObjData *objData,
    struct nestedLitterData *litterData, const std::vector<int> &grpSize
) {
  std::vector<double> Yp = objData->Yp;
  std::vector<double> Yn = objData->Yn;
  std::vector<double> Ls = objData->Ls;
  std::vector<double> Lsc = objData->Lsc;
  std::vector<double> Xi = objData->Xi;
  std::vector<int> Xg = objData->Xg;
  int ngrp = objData->ngrp;
  int Nobs = Yp.size();
  int nparm = parms.size();
  std::vector<double> Ep(Nobs);
  std::vector<double> phi(ngrp);  // phi parms for each group
  std::vector<double> Ysum(Nobs);
  std::vector<double> Ypp(Nobs);
  std::vector<double> Var(Nobs);
  std::vector<double> SR(Nobs);

  // assumes PHI_START = 5
  for (int i = 0; i < ngrp; i++) {
    phi[i] = parms[5 + i];
  }

  // Predict
  Nlogist_Predict(parms, objData, Ep);

  // replace objData with sorted values
  objData->Ls = Ls;
  objData->Lsc = Lsc;
  objData->Yp = Yp;
  objData->Yn = Yn;
  for (int i = 0; i < Nobs; i++) {
    Ysum[i] = Yp[i] + Yn[i];
    Ypp[i] = Ep[i] * Ysum[i];
    Var[i] = Ysum[i] * Ep[i] * ((1 - Ep[i]) * (1.0 + (Ysum[i] - 1.0) * phi[Xg[i]]));
  }

  // compute the goodness of fit test for litter data
  double gfit = 0.0;
  for (int i = 0; i < Nobs; i++) {
    if (Var[i] > 0) {
      SR[i] = Var[i] > 0 ? (Yp[i] - Ypp[i]) / sqrt(Var[i]) : 0;
      gfit += SR[i] * SR[i];
    }
  }
  litterData->chiSq = gfit;

  for (int i = 0; i < Nobs; i++) {
    litterData->dose[i] = Xi[i];
    litterData->LSC[i] = Lsc[i];
    litterData->estProb[i] = Ep[i];
    litterData->litterSize[i] = Ysum[i];
    litterData->expected[i] = Ypp[i];
    litterData->observed[i] = Yp[i];
    litterData->SR[i] = (Var[i] > 0 ? (Yp[i] - Ypp[i]) / sqrt(Var[i]) : 0);
  }
}

void SortNestedData(
    const std::vector<int> &grpSize, std::vector<double> &Xi, std::vector<double> &Ls,
    std::vector<double> &Yp, std::vector<double> &Yn, std::vector<double> &Lsc, bool sortByLsc
) {
  // sorts nested data by Xi -> Ls(or Lsc) -> Yp

  int ngrp = grpSize.size();
  int Nobs = Xi.size();

  std::vector<double> origXi = Xi;
  std::vector<double> origLs = Ls;
  std::vector<double> origYp = Yp;
  std::vector<double> origYn = Yn;
  std::vector<double> origLsc = Lsc;

  std::vector<double> newXi(Nobs);
  std::vector<double> newLs(Nobs);
  std::vector<double> newYp(Nobs);
  std::vector<double> newYn(Nobs);
  std::vector<double> newLsc(Nobs);

  // first pull out each dose group, then sort each dose group, then combine back
  int grpStart = 0;
  for (int i = 0; i < ngrp; i++) {
    int thisGrpSize = grpSize[i];

    std::vector<std::vector<double>> sortV;
    for (int j = grpStart; j < grpStart + thisGrpSize; j++) {
      std::vector<double> tmp;
      tmp.push_back(origXi[j]);
      tmp.push_back(origLs[j]);
      tmp.push_back(origYp[j]);
      tmp.push_back(origYn[j]);
      tmp.push_back(origLsc[j]);
      sortV.push_back(tmp);
    }
    // first sort the grouped vector by the 3rd vector Yp
    // then sort the grouped vector by 2nd vector Ls or 5th vector Lsc
    // then sort the grouped vector by 1st vector Xi
    std::sort(
        sortV.begin(), sortV.end(),
        [](const std::vector<double> &a, const std::vector<double> &b) { return a[2] < b[2]; }
    );
    if (sortByLsc) {
      std::sort(
          sortV.begin(), sortV.end(),
          [](const std::vector<double> &a, const std::vector<double> &b) { return a[4] < b[4]; }
      );
    } else {
      std::sort(
          sortV.begin(), sortV.end(),
          [](const std::vector<double> &a, const std::vector<double> &b) { return a[1] < b[1]; }
      );
    }
    std::sort(sortV.begin(), sortV.end());

    // insert into new vector
    for (int i = grpStart; i < grpStart + thisGrpSize; i++) {
      newXi[i] = sortV[i - grpStart][0];
      newLs[i] = sortV[i - grpStart][1];
      newYp[i] = sortV[i - grpStart][2];
      newYn[i] = sortV[i - grpStart][3];
      newLsc[i] = sortV[i - grpStart][4];
    }

    grpStart += grpSize[i];  // set grpStart index for next group
  }

  // replace original with sorted
  Xi = newXi;
  Ls = newLs;
  Yp = newYp;
  Yn = newYn;
  Lsc = newLsc;
}

void Nlogist_Predict(
    const std::vector<double> &parms, struct nestedObjData *objData, std::vector<double> &P
) {
  std::vector<double> Xi = objData->Xi;
  std::vector<double> Lsc = objData->Lsc;

  double bkg;
  int Nobs = Xi.size();
  for (int i = 0; i < Nobs; i++) {
    bkg = parms[0] + parms[2] * Lsc[i];
    if (Xi[i] <= 0.0) {
      P[i] = bkg;
    } else {
      P[i] =
          bkg + (1.0 - bkg) / (1.0 + exp(-(parms[1] + parms[3] * Lsc[i] + parms[4] * log(Xi[i]))));
    }
  }
}

void Nlogist_vcv(
    std::vector<double> &p, std::vector<bool> &bounded, struct nestedObjData *objData,
    std::vector<std::vector<double>> &vcv
) {
  double tmp;
  int jvar;

  std::vector<bool> Spec = objData->Spec;
  double smin = objData->smin;
  double smax = objData->smax;
  objData->isBMDL = false;

  std::vector<double> ptemp = p;

  if (Spec[0] == Spec[2]) {
    double junk1 = ptemp[0];
    double junk3 = ptemp[2];
    ptemp[0] = junk1 + smin * junk3;
    ptemp[2] = junk1 + smax * junk3;
  }

  // Get a value of h for each parameter
  double hrat = pow(1.0e-16, 0.333333);
  std::vector<double> h(p.size());
  std::vector<double> gradp(p.size());
  std::vector<double> gradm(p.size());

  for (int i = 0; i < ptemp.size(); i++) {
    if (fabs(ptemp[i]) > 1.0e-7) {
      h[i] = hrat * fabs(ptemp[i]);
      tmp = ptemp[i] + h[i];
      h[i] = tmp - ptemp[i];  /// WHY????
    } else {
      h[i] = hrat;
    }
  }

  std::vector<double> saveParms = ptemp;
  int ivar = 0;
  int nvar = 0;

  for (int i = 0; i < ptemp.size(); i++) {
    if (i > 0) saveParms[i - 1] = ptemp[i - 1];
    if (!Spec[i]) {
      nvar++;
      saveParms[i] = ptemp[i] + h[i];
      Nlogist_grad(saveParms, objData, gradp);
      saveParms[i] = ptemp[i] - h[i];
      Nlogist_grad(saveParms, objData, gradm);
      // Now compute the 2nd derivative
      jvar = 0;
      for (int j = 0; j < ptemp.size(); j++) {
        if (!Spec[j]) {
          vcv[ivar][jvar] = -(gradp[jvar] - gradm[jvar]) / (2.0 * h[i]);
          jvar++;
        }
      }
    }
    ivar++;
    nvar++;
  }
}

// Computes the gradient of the nested logistic likelihood function with respect to the user form
// (external of the parameters.  This is to be used in Nlogist_vcv to compute a finite difference
// approximation to the hessian of the likelihood function.
void Nlogist_grad(
    std::vector<double> &p, struct nestedObjData *objData, std::vector<double> &grad
) {
  int nparm = p.size();
  std::vector<double> pint = p;
  // transform the parameters to "internal" form
  for (int j = 5; j < nparm; j++) {
    pint[j] = pint[j] / (1 - pint[j]);
  }

  Nlogist_g(pint, grad, objData);

  // Nlogist_g returns the gradient of the negative loglikelihood
  for (int i = 0; i < nparm; i++) {
    grad[i] *= -1.0;
  }
}

double opt_nlogistic(std::vector<double> &p, struct nestedObjData *objData) {
  std::vector<bool> Spec = objData->Spec;
  int nparm_prior = objData->prior.size() / 2;

  bool fail = true;
  double minf;
  // start with hardcoded parameter limits
  // will move these to header
  double slopeUpperBound = 18.0;
  int nparm = p.size();

  std::vector<double> pbak = p;

  std::vector<double> lb(nparm);
  std::vector<double> ub(nparm);
  // alpha
  lb[0] = objData->prior[0];            // 0.0;
  ub[0] = objData->prior[nparm_prior];  // 1.0;
  // beta
  lb[1] = objData->prior[1];                //-1.0*DBL_MAX;
  ub[1] = objData->prior[nparm_prior + 1];  // DBL_MAX;
  // theta1
  lb[2] = objData->prior[2];  // 0.0;
  if (Spec[2]) {
    ub[2] = 0.0;
  } else {
    ub[2] = objData->prior[nparm_prior + 2];  // 1.0;
  }
  // theta2
  if (Spec[3]) {
    lb[3] = 0.0;
    ub[3] = 0.0;
  } else {
    lb[3] = objData->prior[3];                //-1.0*DBL_MAX;
    ub[3] = objData->prior[nparm_prior + 3];  // DBL_MAX;
  }
  // rho
  lb[4] = objData->prior[4];                // 0.0/1.0 (unrestricted/restricted)
  ub[4] = objData->prior[nparm_prior + 4];  // slopeUpperBound;
  // phi(s)
  for (int i = 5; i < nparm; i++) {
    lb[i] = objData->prior[5];  // 0.0;
    if (Spec[i]) {
      ub[i] = 0.0;
    } else {
      ub[i] = objData->prior[nparm_prior + 5];  // DBL_MAX;
    }
  }

  int result;

  nlopt::opt opt(nlopt::LN_AUGLAG, nparm);

  nlopt::opt local_opt(nlopt::LD_LBFGS, nparm);
  nlopt::opt local_opt2(nlopt::LN_SBPLX, nparm);

  local_opt.set_xtol_abs(1e-3);
  local_opt.set_initial_step(1e-4);
  local_opt.set_maxeval(10000);

  local_opt2.set_xtol_abs(1e-3);
  local_opt2.set_initial_step(1e-4);
  local_opt2.set_maxeval(10000);

  local_opt.set_lower_bounds(lb);
  local_opt2.set_lower_bounds(lb);
  local_opt.set_upper_bounds(ub);
  local_opt2.set_upper_bounds(ub);

  // Description of optimization problem
  // minimize -log-likelihood
  // inequality constraints
  //	alpha + theta1*rij >=0
  //	alpha + theta1*rij < 1

  if (Spec[0] == Spec[3]) {
    // set inequality constraint alpha + Theta1*Sij>=0
    // inequality restraint not compatible with LD_LBFGS
    // opt= nlopt::opt(nlopt::LD_SLSQP, nparm);
    opt.add_inequality_constraint(nestedInequalityConstraint, &objData, 1e-8);
  }
  bool good_opt = false;
  int opt_iter = 0;

  std::vector<double> init(nparm);
  for (int i = 0; i < nparm; i++) init[i] = 1e-4;
  opt.set_initial_step(init);
  local_opt.set_initial_step(init);

  opt.set_min_objective(objfunc_nlogistic_ll, objData);
  while (opt_iter < 2 && !good_opt) {
    opt_iter++;
    if (opt_iter == 0)
      opt.set_local_optimizer((const nlopt::opt)local_opt);
    else
      opt.set_local_optimizer((const nlopt::opt)local_opt2);

    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);
    opt.set_xtol_abs(1e-4);
    opt.set_maxeval(20000);
    // Ensure that starting values are within bounds
    for (int i = 0; i < nparm; i++) {
      double temp = p[i];
      if (temp < lb[i])
        temp = lb[i];
      else if (temp > ub[i])
        temp = ub[i];
      p[i] = temp;
    }  // end for

    try {
      result = opt.optimize(p, minf);
      good_opt = true;  // optimization succeded
    } catch (nlopt::roundoff_limited &exec) {
      good_opt = false;
      std::cout << "Round Off Limited" << std::endl;
    } catch (nlopt::forced_stop &exec) {
      good_opt = false;
      std::cout << "Forced Stop" << std::endl;
    } catch (const std::invalid_argument &exc) {
      std::cout << "Invalid Argument" << std::endl;
      good_opt = false;
    } catch (const std::exception &exc) {
      std::cout << "Std Exception" << std::endl;
      good_opt = false;
    } catch (...) {
      std::cout << "Default Exception" << std::endl;
      good_opt = false;
    }

    if (result > 5) {  // Either 5 =
      good_opt = false;
    }
  }

  objData->xlk = -1.0 * minf;

  return result;
}

void probability_inrange(double *ex) {
  if (*ex < 1.0e-7) *ex = 1.0e-7;
  if (*ex > 0.9999999) *ex = 0.9999999;
}

double objfunc_nlogistic_ll(const std::vector<double> &p, std::vector<double> &grad, void *data) {
  nestedObjData *objData = reinterpret_cast<nestedObjData *>(data);
  // from data struct
  std::vector<double> Ls = objData->Ls;
  std::vector<double> Xi = objData->Xi;
  std::vector<int> Xg = objData->Xg;
  std::vector<double> Yp = objData->Yp;
  std::vector<double> Yn = objData->Yn;
  double smax = objData->smax;
  double smin = objData->smin;
  bool isBMDL = objData->isBMDL;
  double D = objData->tD;
  double sijfixed = objData->sijfixed;
  int riskType = objData->riskType;
  double BMR = objData->BMR;

  if (!grad.empty()) {
    Nlogist_g(p, grad, objData);
  }

  // negative log-likelihood calc
  double loglike = Nlogist_lk(p, objData);

  return loglike;
}

double nestedInequalityConstraint(
    const std::vector<double> &x, std::vector<double> &grad, void *data
) {
  // alpha + theta1 *Rij > = 0
  nestedObjData *objData = reinterpret_cast<nestedObjData *>(data);
  double sijfixed = objData->sijfixed;

  if (!grad.empty()) {
    for (size_t i = 0; i < grad.size(); i++) {
      grad[i] = 0.0;
    }
    grad[0] = -1.0 * x[2];
    grad[2] = -1.0 * x[0];
  }

  double constraint = -1.0 * (x[0] + x[2] * sijfixed);
  return constraint;
}

// Used to comput the log-likelihood for Nlogistic model
double Nlogist_lk(std::vector<double> p, struct nestedObjData *objData) {
  std::vector<double> Yp = objData->Yp;
  std::vector<double> Yn = objData->Yn;
  std::vector<int> Xg = objData->Xg;

  int Nobs = Yp.size();
  int nparms = p.size();
  std::vector<double> probs(Nobs);
  std::vector<std::vector<double>> gradij(Nobs, std::vector<double>(5));

  bool compgrad = false;
  Nlogist_probs(probs, p, compgrad, gradij, objData);

  double tm, tm1, tm2, tm3;
  int plus4, j;
  double xlk = 0.0;
  for (int i = 0; i < Nobs; i++) {
    tm1 = 0.0;
    tm2 = 0.0;
    tm3 = 0.0;
    plus4 = 4 + Xg[i];
    j = (int)Yp[i];
    if (probs[i] < CloseToZero && j > 0) {
      tm1 -= 40.0;
    } else {
      for (int k = 1; k <= j; k++) {
        tm = probs[i] + (k - 1) * p[plus4];
        tm1 += log(tm);
      }
    }
    j = (int)Yn[i];
    if (probs[i] >= 1.0 && j > 0) {
      tm2 -= 40.0;
    } else {
      for (int k = 1; k <= j; k++) {
        tm = 1.0 - probs[i] + (k - 1) * p[plus4];
        tm2 += log(tm);
      }
    }
    j = (int)(Yn[i] + Yp[i]);
    for (int k = 1; k <= j; k++) {
      tm = 1.0 + (k - 1) * p[plus4];
      tm3 += log(tm);
    }
    xlk += (tm1 + tm2 - tm3);
  }

  // returns negative log likelihood
  return -1.0 * xlk;
}

// Used to compute the gradients for Nlogist_model.  Parameters are in "internal" transformed form.
double Nlogist_g(std::vector<double> p, std::vector<double> &g, struct nestedObjData *objData) {
  std::vector<double> Yp = objData->Yp;
  std::vector<double> Yn = objData->Yn;
  std::vector<int> Xg = objData->Xg;
  std::vector<bool> Spec = objData->Spec;

  int Nobs = Yp.size();
  int nparm = p.size();
  std::vector<double> probs(Nobs);
  std::vector<std::vector<double>> gradij(Nobs, std::vector<double>(5));
  double ex, tm, tm1, tm2, tm3, tm1a, tm2a, tm3a, tm12;
  int plus5, j;
  std::vector<double> tmp_g(nparm);
  std::vector<double> dd(nparm);

  bool compgrad = true;
  Nlogist_probs(probs, p, compgrad, gradij, objData);

  // initial tmp g[j]'s
  for (int i = 0; i < nparm; i++) {
    tmp_g[i] = 0.0;
  }
  for (int i = 0; i < Nobs; i++) {
    ex = probs[i];

    // compute first partial derivatives
    tm1 = tm2 = tm3 = 0.0;
    tm1a = tm2a = tm3a = 0.0;
    for (int j = 5; j < nparm; j++) {
      dd[j] = 0.0;
    }
    plus5 = 4 + Xg[i];
    j = (int)Yp[i];
    if (ex > 0.0) {
      for (int k = 1; k <= j; k++) {
        tm = ex + (k - 1) * p[plus5];
        tm1 += 1.0 / tm;
        tm1a += (1.0 / tm) * (k - 1);
      }
    }
    j = (int)Yn[i];
    if (ex < 1.0) {
      for (int k = 1; k <= j; k++) {
        tm = 1.0 - ex + (k - 1) * p[plus5];
        tm2 += 1.0 / tm;
        tm2a += (1.0 / tm) * (k - 1);
      }
    }
    j = (int)(Yn[i] + Yp[i]);
    for (int k = 1; k <= j; k++) {
      tm = 1.0 + (k - 1) * p[plus5];
      if (tm == 0.0) tm = 0.000001;
      tm3 += 1.0 / tm;
      tm3a += (1.0 / tm) * (k - 1);
    }
    tm12 = (tm1 - tm2);
    for (int j = 0; j < 5; j++) {
      dd[j] = gradij[i][j] * tm12;
    }
    dd[plus5] = (tm1a + tm2a - tm3a);
    for (int j = 0; j < nparm; j++) {
      tmp_g[j] -= dd[j];
    }
  }
  // end of 1st partial derivative
  g = tmp_g;
  for (int j = 0; j < nparm; j++) {
    if (Spec[j]) {
      g[j] = 0.0;
    }
  }

  return 0;
}

void Nlogist_probs(
    std::vector<double> &probs, const std::vector<double> &p, bool compgrad,
    std::vector<std::vector<double>> &gradij, struct nestedObjData *objData
) {
  bool isBMDL = objData->isBMDL;
  double smax = objData->smax;
  double smin = objData->smin;
  std::vector<double> Ls = objData->Ls;
  std::vector<double> Xi = objData->Xi;
  double sijfixed = objData->sijfixed;
  int riskType = objData->riskType;
  double BMR = objData->BMR;
  double tD = objData->tD;
  std::vector<bool> Spec = objData->Spec;

  double spij, smij, ex, ex1, ex2, ex3, ex4, dd2;

  double sdif = smax - smin;
  double spfixed = (smax - sijfixed) / sdif;
  double snfixed = (sijfixed - smin) / sdif;
  int nparms = p.size();
  int Nobs = Xi.size();

  std::vector<double> pint(nparms);
  for (int i = 0; i < nparms; i++) {
    pint[i] = p[i];
  }

  if (isBMDL) {
    if (Spec[0] == Spec[2]) {
      if (riskType == 1) {  // Extra
        pint[1] = log(BMR / (1.0 - BMR)) - pint[3] * sijfixed - pint[4] * log(tD);
      } else {  // Added
        pint[1] = -log((1 - pint[0] * spfixed - pint[2] * snfixed) / BMR - 1) - pint[3] * sijfixed -
                  pint[4] * log(tD);
      }
    } else {
      if (riskType == 1) {  // Extra
        pint[1] = log(BMR / (1.0 - BMR)) - pint[3] * sijfixed - pint[4] * log(tD);
      } else {  // Added
        pint[1] = -log((1 - pint[0] - pint[2] * sijfixed) / BMR - 1) - pint[3] * sijfixed -
                  pint[4] * log(tD);
      }
    }
  }

  for (int i = 0; i < Nobs; i++) {
    spij = (smax - Ls[i]) / sdif;
    smij = (Ls[i] - smin) / sdif;

    // enforce alpha+Theta1*Sij >= 0
    if (Spec[0] == Spec[2]) {
      ex = spij * pint[0] + smij * pint[2];
    } else {
      ex = pint[0] + Ls[i] * pint[2];
    }
    ex2 = 1.0 - ex;
    ex1 = 0.0;
    if (Xi[i] > 0) {
      ex1 = exp(-(pint[1] + pint[3] * Ls[i] + pint[4] * log(Xi[i])));
      ex = ex + ex2 / (1 + ex1);
    }
    probability_inrange(&ex);
    probs[i] = ex;

    if (compgrad) {
      ex3 = ex1 / ((1.0 + ex1) * (1.0 + ex1));
      dd2 = ex2 * ex3;

      if (isBMDL) {
        if (riskType == 1) {  // Extra
          ex4 = 0.0;
        } else {  // Added
          ex4 = 1.0 / (ex2 - BMR);
        }
        if (Spec[0] == Spec[2]) {
          if (Xi[i] > 0) {
            gradij[i][0] = (spij * (1.0 - 1.0 / (1.0 + ex1)) + dd2 * ex4 * spfixed);
            gradij[i][2] = (smij * (1.0 - 1.0 / (1.0 + ex1)) + dd2 * ex4 * snfixed);
          } else {
            gradij[i][0] = spij;
            gradij[i][2] = smij;
          }
        } else {
          if (Xi[i] > 0) {
            gradij[i][0] = 1.0 - 1.0 / (1.0 + ex1) + dd2 * ex4;
            gradij[i][2] = Ls[i] * gradij[i][0];
          } else {
            gradij[i][0] = 1.0;
            gradij[i][2] = Ls[i];
          }
        }
      } else {
        // enforce alpha+Theta1*Sij >= 0
        if (Spec[0] == Spec[2]) {
          if (Xi[i] > 0) {
            gradij[i][0] = spij * (1.0 - 1.0 / (1.0 + ex1));
            gradij[i][2] = smij * (1.0 - 1.0 / (1.0 + ex1));
          } else {
            // Case where Xi[i] == 0
            gradij[i][0] = spij;
            gradij[i][2] = smij;
          }
        } else {
          if (Xi[i] > 0) {
            gradij[i][0] = 1.0 - 1.0 / (1.0 + ex1);
            gradij[i][2] = Ls[i] * gradij[i][0];
          } else {
            // Case where Xi[i] == 0
            gradij[i][0] = 1.0;
            gradij[i][2] = Ls[i];
          }
        }  // end of derivates that depend on the transformation of alpha and theta1
        if (Xi[i] > 0) {
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
}

void Nlogist_BMD(
    struct python_nested_analysis *pyAnal, struct python_nested_result *pyRes, double smin,
    double smax, double sijfixed, double xmax, struct nestedObjData *objData
) {
  std::vector<double> p = pyRes->parms;
  std::vector<bool> Spec = objData->Spec;
  int nparm = p.size();
  double BMR = pyAnal->BMR;
  int riskType = pyAnal->BMD_type;
  double junk1, junk3, ck;
  double CL = 1.0 - pyAnal->alpha;

  // If ML is the value of the maximized log-likelihood, then ML - LR is the value log-likehood at
  // the BMDL or BMDU
  if (CL < 0.5) {
    objData->LR = QCHISQ(1.0 - 2 * CL, 1) / 2.0;
  } else {
    objData->LR = QCHISQ(2 * CL - 1.0, 1) / 2.0;
  }

  std::vector<double> pint = p;

  // transform parameters into "internal" form
  for (int i = 5; i < nparm; i++) {
    pint[i] = pint[i] / (1 - pint[i]);  // Phi->Psi
  }

  if (Spec[0] == Spec[2]) {
    junk1 = pint[0];
    junk3 = pint[2];
    pint[0] = junk1 + smin * junk3;
    pint[2] = junk1 + smax * junk3;
  }

  double sdif = smax - smin;

  double spfixed = (smax - sijfixed) / sdif;
  double snfixed = (sijfixed - smin) / sdif;

  if (riskType == 1) {
    // extra risk
    ck = -log((1 - BMR) / BMR);
  } else {
    // added risk
    if (Spec[0] == Spec[2]) {
      ck = -1.0 * log((1 - pint[0] * spfixed - pint[2] * snfixed) / BMR - 1);
    } else {
      ck = -1.0 * log((1 - pint[0] - pint[2] * sijfixed) / BMR - 1);
    }
  }
  objData->ck = ck;

  double BMD;
  if (pint[4] <= (ck - pint[1] - pint[3] * sijfixed) / 250) {
    // Power parameter is essentially zero.  BMD is set to 100  * max(Dose)
    BMD = 100 * xmax;
  } else {
    BMD = exp((ck - pint[1] - pint[3] * sijfixed) / pint[4]);
  }

  pyRes->bmdsRes.BMD = BMD;
  pyRes->bmd = BMD;

  // Search for BMDL
  objData->isBMDL = true;
  double stepsize = 0.5;  // Start close to the BMD and work downwards
  double xb = BMD;
  double xa = xb * stepsize;
  double tol = max(BMD * 0.001, 0.0000001);
  double fb = DBL_MAX;

  pyRes->bmdsRes.BMDL = calcNlogisticCLs(xa, xb, pint, objData, true);

  // Search for BMDU

  objData->isBMDL = true;
  stepsize = 0.5;  // Start close to the BMD and work upwards
  xb = BMD;
  xa = xb * (1.0 + stepsize);
  tol = max(BMD * 0.001, 0.0000001);
  fb = DBL_MAX;

  pyRes->bmdsRes.BMDU = calcNlogisticCLs(xa, xb, pint, objData, false);
}

double calcNlogisticCLs(
    double xa, double xb, std::vector<double> &pint, struct nestedObjData *objData, bool isLower
) {
  double stepsize = 0.5;  // Start close to the BMD and work upwards
  double tol = max(xb * 0.001, 0.0000001);
  double fb = DBL_MAX;

  int nparm = pint.size();
  std::vector<double> pa(nparm);
  std::vector<double> pb(nparm);

  for (int i = 0; i < nparm; i++) {
    pa[i] = pb[i] = pint[i];
  }

  double fa = BMDL_func(pa, xa, tol, objData);
  if (!isLower) fa *= -1;

  // Look for a value of xa on the other side of the BMDL.  We know we're there when fa > 0.
  // Stop if xa gets too small, or the profile likelihood gets flat (fabs(fa - fb) too small).

  int trip = 0;
  double tmp = fabs(fa - fb);
  while (fa < 0.0 && xa > DBL_MIN && fabs(fa - fb) > DBL_EPSILON) {
    xb = xa;
    fb = fa;
    for (int i = 0; i < nparm; i++) pb[i] = pa[i];
    xa *= stepsize;

    fa = BMDL_func(pa, xa, tol, objData);
    if (!isLower) fa *= -1;
    trip++;
  }

  double val;  // either BMDL or BMDU
  if (fa < 0.0) {
    val = BMDS_MISSING;
  } else {
    val = zeroin_nested(xa, xb, 1.0e-10, BMDL_func, pb, 1.0e-14, objData);
  }

  return val;
}

// QCHISQ - inverse chi-square function
double QCHISQ(double p, int m) {
  double df = (double)m;
  double x = gsl_cdf_chisq_Pinv(p, df);
  return x;
}

double CHISQ(double x, int m) {
  double df = (double)m;
  double p = gsl_cdf_chisq_P(x, df);
  return p;
}

void outputObjData(struct nestedObjData *objData) {
  std::cout << "nestedObjData-------------" << std::endl;
  std::cout << "ck:" << objData->ck << std::endl;
  std::cout << "LR:" << objData->LR << std::endl;
  std::cout << "xlk:" << objData->xlk << std::endl;
  std::cout << "BMD_lk:" << objData->BMD_lk << std::endl;
  std::cout << "tD:" << objData->tD << std::endl;
  std::cout << "sijfixed:" << objData->sijfixed << std::endl;
  std::cout << "riskType:" << objData->riskType << std::endl;
  std::cout << "BMR:" << objData->BMR << std::endl;

  for (int i = 0; i < objData->Spec.size(); i++) {
    std::cout << "i:" << i << ", Spec:" << objData->Spec[i] << std::endl;
  }
}

// BMDL_func - used to compare the values of functions BMDL_f (the X^2 value) at the point D,
//    given the parm p[] and the number of parm.  Input parameters are in the "internal" form.
//    This routine is called by zeroin()
double BMDL_func(std::vector<double> &p, double D, double gtol, struct nestedObjData *objData) {
  double fD;
  int junk;

  objData->tD = D;
  objData->tol = gtol;
  // std::vector<double> parms = p;
  double retVal = opt_nlogistic(p, objData);

  // set result parms to return array
  // p = parms;
  fD = objData->BMD_lk - objData->xlk - objData->LR;
  return fD;
}

double round_to(double value, double precision) {
  double res = std::round(value / precision) * precision;
  return res;
}

/*
 *  ************************************************************************
 *	    		    C math library
 * function ZEROIN_NESTED - obtain a function zero within the given range
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

double zeroin_nested(
    double ax, double bx, double tol,
    double (*f)(std::vector<double> &, double, double, struct nestedObjData *),
    std::vector<double> &Parms, double ck, struct nestedObjData *objData
)
/* ax        Left border | of the range */
/* bx        Right border | the root is sought*/
/* f	  Function under investigation */
/* nparm     number of parameters in Parms */
/* Parms     vector of parameters to pass to f */
/* tol       Acceptable tolerance for the root */
/* gtol      tolerance to pass to f */
// TODO: clean up comments after debugging is complete
{
  double a, b, c; /* Abscissae, descr. see above	*/
  double fa;      /* f(a)				*/
  double fb;      /* f(b)				*/
  double fc;      /* f(c)				*/

  int nparm = Parms.size();

  a = ax;
  b = bx;
  fa = (*f)(Parms, a, ck, objData);
  fb = (*f)(Parms, b, ck, objData);
  c = a;
  fc = fa;
  int pass = 1;
  int maxPass = 10000;
  for (;;) /* Main iteration loop	*/
  {
    double prev_step = b - a; /* Distance from the last but one*/
                              /* to the last approximation	*/
    double tol_act;           /* Actual tolerance		*/
    double p;                 /* Interpolation step is calcu- */
    double q;                 /* lated in the form p/q; divi- */
                              /* sion operations is delayed   */
    /* until the last moment	*/
    double new_step;           /* Step at this iteration       */
    if (fabs(fc) < fabs(fb)) { /* Swap data for b to be the 	*/
      a = b;
      b = c;
      c = a; /* best approximation		*/
      fa = fb;
      fb = fc;
      fc = fa;
    }
    tol_act = 2 * DBL_EPSILON * fabs(b) + tol / 2;
    new_step = (c - b) / 2;

    if (fabs(new_step) <= tol_act || fb == (double)0) {
      return b; /* Acceptable approx. is found	*/
    }

    /* Decide if the interpolation can be tried	*/
    if (fabs(prev_step) >= tol_act /* If prev_step was large enough*/
        && fabs(fa) > fabs(fb))    /* and was in true direction,	*/
    {                              /* Interpolatiom may be tried	*/
      double t1, cb, t2;
      cb = c - b;
      if (a == c)     /* If we have only two distinct	*/
      {               /* points linear interpolation 	*/
        t1 = fb / fa; /* can only be applied		*/
        p = cb * t1;
        q = 1.0 - t1;
      } else /* Quadric inverse interpolation*/
      {
        q = fa / fc;
        t1 = fb / fc;
        t2 = fb / fa;
        p = t2 * (cb * q * (q - t1) - (b - a) * (t1 - 1.0));
        q = (q - 1.0) * (t1 - 1.0) * (t2 - 1.0);
      }
      if (p > (double)0) /* p was calculated with the op-*/
        q = -q;          /* posite sign; make p positive	*/
      else               /* and assign possible minus to	*/
        p = -p;          /* q				*/

      if (p < (0.75 * cb * q - fabs(tol_act * q) / 2) /* If b+p/q falls in [b,c]*/
          && p < fabs(prev_step * q / 2))             /* and isn't too large	*/
        new_step = p / q;                             /* it is accepted	*/
                                                      /* If p/q is too large then the	*/
                                                      /* bissection procedure can 	*/
                                                      /* reduce [b,c] range to more	*/
                                                      /* extent			*/
    }

    if (fabs(new_step) < tol_act) /* Adjust the step to be not less*/
    {
      if (new_step > (double)0) /* than tolerance		*/
        new_step = tol_act;
      else
        new_step = -tol_act;
    }

    a = b;
    fa = fb; /* Save the previous approx.	*/
    b += new_step;
    fb = (*f)(Parms, b, ck, objData);               /* Do step to a new approxim.	*/
    if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) { /* Adjust c for it to have a sign*/
      c = a;
      fc = fa; /* opposite to that of b	*/
    }
    pass++;
    if (pass > maxPass) return BMDS_MISSING;
  }
}
