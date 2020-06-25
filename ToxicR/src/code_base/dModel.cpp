#include "stdafx.h"
#ifdef R_COMPILATION
    //necessary things to run in R    
    #include <RcppEigen.h>
    #include <RcppGSL.h>
#else 
    #include <Eigen/Dense>
#endif

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include "bmd_calculate.h"
#include "dModel.h"
#include <iostream>



using namespace std;
using namespace Eigen;

dModel::dModel() {}

dModel::~dModel() {}

void dModel::calcGoF(dGoF_t * zOut, BMDSInputData_t * dataIn, const int nRows, int npEst)
{
  double chisq = 0;
  int df = nRows - npEst;

  for (int i = 0; i < nRows; i++) {
    GoFRow_t *pzRow = &zOut->pzRow[i];
    double residual, stdErr;
    double scaledResidual = 0.0;
    double N = dataIn[i].groupSize; // # subjects in dose group
    double Yp = dataIn[i].response; // Observed # postive responses
    double Yn = N - Yp; // Observed # negative responses

    pzRow->dose = dataIn[i].dose;
    pzRow->estProb = this->mean(dataIn[i].dose);
    pzRow->observed = Yp;
    pzRow->size = N;
    pzRow->expected = pzRow->estProb * N;
    residual = pzRow->observed - pzRow->expected;
    stdErr = sqrt(pzRow->expected * (1.0 - pzRow->estProb));
    if (stdErr > 0) scaledResidual = residual / stdErr;
    chisq += scaledResidual * scaledResidual;
    pzRow->scaledResidual = scaledResidual;

    // Calculate error bar endpoints using Fleiss confidence limits
    double pHat = Yp / N; // Observed probability
    constexpr auto myAlpha = 0.05; // Alpha value for 95% confidence limit
    double z = gsl_cdf_ugaussian_Pinv(1.0 - myAlpha / 2); // Z score
    pzRow->ebLower = (2 * Yp + z * z - 1) - z * sqrt(z*z - (2 + 1 / N) + 4 * pHat*((N - Yp) + 1));
    pzRow->ebLower /= 2.0 * (N + z * z);
    pzRow->ebUpper = (2 * Yp + z * z + 1) + z * sqrt(z*z + (2 - 1 / N) + 4 * pHat*((N - Yp) - 1));
    pzRow->ebUpper /= 2.0 * (N + z * z);

  } // for
  zOut->chiSquare = chisq;
  zOut->df = df;
  zOut->pvalue = 1.0 - gsl_cdf_chisq_P(chisq, df);
  zOut->n = nRows;

  return;
}

bmds_loglogistic::bmds_loglogistic() {}

bmds_loglogistic::~bmds_loglogistic() {}

void bmds_loglogistic::setParms(double * parms)
{
  g = parms[0];
  a = parms[1];
  b = parms[2];
}

double bmds_loglogistic::mean(double dose)
{
  double result = g;
  if (dose <= 0.0) goto bailout;
  result = g + (1 - g) / (1 + exp(-a - b * log(dose)));

bailout:
  return result;
}


bmds_gamma::bmds_gamma() {}

bmds_gamma::~bmds_gamma() {}

void bmds_gamma::setParms(double * parms)
{
  g = parms[0];
  a = parms[1];
  b = parms[2];
}

double bmds_gamma::mean(double dose)
{
  double result = g;
  if (dose <= 0.0) goto bailout;
  result = g + (1 - g)*(gsl_cdf_gamma_P(b*dose, a, 1));

bailout:
  return result;
}

bmds_dhill::bmds_dhill() {}

bmds_dhill::~bmds_dhill() {}

void bmds_dhill::setParms(double * parms)
{
  g = parms[0];
  n = parms[1];
  a = parms[2];
  b = parms[3];
}

double bmds_dhill::mean(double dose)
{
  double result = g;
  if (dose <= 0.0) goto bailout;
  result = g + (1 - g)*n / (1 + exp(-a - b * log(dose)));

bailout:
  return result;
}

bmds_logistic::bmds_logistic() {}

bmds_logistic::~bmds_logistic() {}

void bmds_logistic::setParms(double * parms)
{
  a = parms[0];
  b = parms[1];
}

double bmds_logistic::mean(double dose)
{
  double result;
  result = 1 / (1 + exp(-a - b * dose));
  return result;
}

bmds_logprobit::bmds_logprobit() {}

bmds_logprobit::~bmds_logprobit() {}

void bmds_logprobit::setParms(double * parms)
{
  g = parms[0];
  a = parms[1];
  b = parms[2];
}

double bmds_logprobit::mean(double dose)
{
  double result = g;
  if (dose <= 0.0) goto bailout;
  result = g + (1 - g)*gsl_cdf_gaussian_P(a + b * log(dose), 1);

bailout:
  return result;
}

bmds_probit::bmds_probit() {}

bmds_probit::~bmds_probit() {}

void bmds_probit::setParms(double * parms)
{
  a = parms[0];
  b = parms[1];
}

double bmds_probit::mean(double dose)
{
  double result;
  result = gsl_cdf_gaussian_P(a + b * dose, 1.0);
  return result;
}

bmds_qlinear::bmds_qlinear() {}

bmds_qlinear::~bmds_qlinear() {}

void bmds_qlinear::setParms(double * parms)
{
  g = parms[0];
  b = parms[1];
}

double bmds_qlinear::mean(double dose)
{
  double result = g;
  if (dose <= 0.0) goto bailout;
  result = g + (1.0 - g)*(1.0 - exp(-b * dose));

bailout:
  return result;
}

bmds_weibull::bmds_weibull() {}

bmds_weibull::~bmds_weibull() {}

void bmds_weibull::setParms(double * parms)
{
  g = parms[0];
  a = parms[1];
  b = parms[2];
}

double bmds_weibull::mean(double dose)
{
  double result = g;
  if (dose <= 0.0) goto bailout;
  result = g + (1 - g) * (1 - exp(-b * pow(dose, a)));

bailout:
  return result;
}

bmds_multistage::bmds_multistage(int val)
{
  degree = val;
  nParms = val + 1;
  p = new double[nParms];
}

bmds_multistage::~bmds_multistage()
{
  delete[] p;
}

void bmds_multistage::setParms(double * parms)
{
  for (int i = 0; i < nParms; i++) {
    p[i] = parms[i];
  }
}

double bmds_multistage::mean(double dose)
{
  double result = p[nParms - 1] * dose;
  for (int i = nParms - 2; i > 0; i--) {
    result = (result + p[i]) * dose;
  }
  result = p[0] + (1.0 - p[0]) * (1.0 - exp(-result));
  return result;
}
