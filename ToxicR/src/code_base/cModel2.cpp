#include "stdafx.h"

#ifdef R_COMPILATION
    //necessary things to run in R    
    #include <RcppEigen.h>
    #include <RcppGSL.h>
#else 
    #include <Eigen/Dense>
#endif

#include "bmds_entry.h"
#include "bmds_dmodels.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_cdf.h"
//#include "bmds_cmodel_entry.h"
#include <vector>
#include "cModel.h"
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;
using namespace Eigen;

#define LOGFILENAME __FUNCTION__ ## ".log"

cModel::cModel() {}

cModel::~cModel() {}

void cModel::setPriors(Eigen::MatrixXd priors)
{
  zpriors = priors;
} // cModel::setPriors

int cModel::calcGoF(cGoFRow_t * zOut, BMDSInputType_t inputType, BMDSInputData_t * dataIn,
  const int nRows, int npEst, Eigen::MatrixXd X, Eigen::MatrixXd Y)
{
  //cout << "calGoF: bCV= " << bConstVar << " rho= " << rho << " lnalpha= " << lnalpha << endl;
  int nGroups = 0; // Number of unique dose groups
  BMDSInputData_t *pData = NULL;

#if defined(_DEBUG) || defined(DEBUGLOG)
  ofstream file;
  file.open("calcGoF.log" /*, fstream::app */);
  file << "Data inputType= " << inputType << endl;
#endif
  // Rollup the dosegroups if individual responses given
  // ** ASSUME that rows are ordered by dose!!!
  if (inputType == eCont_2) {

    // For simplicity allocate nRows of space to store the summarized values
    // since there usually aren't many rows.
    pData = new BMDSInputData_t[nRows];
    double gSize = 1;
    double gDose = dataIn[0].dose;
    double gSum = dataIn[0].response; // Sum of responses for the current group
    double gSum2 = dataIn[0].response * dataIn[0].response; // sum of squared response
    int g = 0; // Index of current dose group
#if defined(_DEBUG) || defined(DEBUGLOG)
    file << "Summarizing individual response data...\n" << endl;
    file << "#   Dose\tSize\tMean\tALT StDEV" << endl;
#endif
    for (int i = 1; i < nRows; i++) { // start at 2nd element
      if (gDose != dataIn[i].dose) { // New dose group
        pData[g].dose = gDose;
        pData[g].groupSize = gSize;
        pData[g].response = gSum / gSize;
        pData[g].col4 = sqrt((gSum2 - gSum * gSum / gSize) / (gSize - 1));
#if defined(_DEBUG) || defined(DEBUGLOG)
        file << g << "  " << gDose << '\t' << gSize << '\t'
          << pData[g].response << '\t' << pData[g].col4 << '\t' << endl;
#endif
        // reset variables for next group
        g++;
        gSize = 1;
        gDose = dataIn[i].dose;
        gSum = dataIn[i].response;
        gSum2 = dataIn[i].response * dataIn[i].response;
      }
      else {
        gSize++;
        gSum += dataIn[i].response;
        gSum2 += dataIn[i].response * dataIn[i].response;
      } // end if (gDose != dataIn[i].dose)
    } // for
    // Process the last dose group
    pData[g].dose = gDose;
    pData[g].groupSize = gSize;
    pData[g].response = gSum / gSize;
    pData[g].col4 = sqrt((gSum2 - gSum * gSum / gSize) / (gSize - 1));
#if defined(_DEBUG) || defined(DEBUGLOG)
    file << g << "  " << gDose << '\t' << gSize << '\t'
      << pData[g].response << '\t' << pData[g].col4 << '\t' << endl << endl;
#endif
    nGroups = g + 1; // +1 b/c 'g' is a zero-based index

    //// Do a 2nd pass to calculate the st. dev
    //int j = 0; // index for individual data rows
    //for (g = 0; g < nGroups; g++) {
    //  double mean = pData[g].response;
    //  double N = pData[g].groupSize;
    //  double s2 = 0; // temp variance variable
    //  for (int i = 0; i < N; i++) {
    //    s2 += (dataIn[j].response - mean)*(dataIn[j].response - mean)
    //      / (N-1);
    //    j++;
    //  } // for i
    //  pData[g].col4 = sqrt(s2);
    //} // for g
#if defined(_DEBUG) || defined(DEBUGLOG)
    file << "Numer of unique dose groups= " << nGroups << endl << endl;
#endif
  }
  else {
    pData = dataIn;
    nGroups = nRows;
  } // end if (inputType == eCont_2)

  //int df = nGroups - npEst;
#if defined(_DEBUG) || defined(DEBUGLOG)
  file << "\nlnalpha= " << lnalpha << " rho= " << rho << endl;
  file << "Dose\tSize\tMean\tStDev\tEstMean\tEstStd\tScRes\tebLower\tebUpper" << endl;
#endif

  // Always use observed variance for error bars
  //// If constant variance, calculate common variance for error bars
  //double sp = 0.0, nTot = 0.0;
  //if (bConstVar) {
  //  for (int i = 0; i < nGroups; i++) {
  //    double size = pData[i].groupSize;
  //    nTot += size;
  //    sp += (size - 1) * pData[i].col4 * pData[i].col4;
  //  } // for g
  //  sp /= (nTot - nGroups);
  //} // if (bConstVar)

  /**********     Calculate GoF values plus error bars     **********/

  for (int i = 0; i < nGroups; i++) {
    cGoFRow_t *pzRow = &zOut[i];
    double size = pData[i].groupSize;
    double estMean, estStDev;
    // For lognormal, calculate the median and geometric st. dev. for orig data
    double calcMedian = BMDS_BLANK_VALUE, calcGSD = BMDS_BLANK_VALUE;
    if (bLognormal) {
      estMean = this->mean(pData[i].dose);
      estStDev = exp(sqrt(exp(lnalpha + rho * log(abs(estMean)))));
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // These values will be slightly off for individual d-r data, but close
      // enough fow now. This should be reworked to handle individual data.
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double x = pData[i].response;
      double y = pData[i].col4;
      calcMedian = exp(log(x) - log(1 + pow(y / x, 2.0)) / 2);
      calcGSD = exp(sqrt(log(1.0 + pow(y / x, 2.0))));
    }
    else {
      estMean = this->mean(pData[i].dose);
      //double estStDev = sqrt(alpha * pow(abs(estMean), rho));
      estStDev = sqrt(exp(lnalpha + rho * log(abs(estMean))));
      calcMedian = pData[i].response;
      calcGSD = pData[i].col4;
    } // end if (bLognormal)
    double scaledResidual = 0.0;
    const double nearZero = 0.00000001; // Value used for consistency with BMDS 2.7

    pzRow->dose = pData[i].dose;
    pzRow->size = size;
    pzRow->obsMean = pData[i].response;
    pzRow->obsStDev = pData[i].col4;
    pzRow->estMean = estMean;
    pzRow->estStDev = estStDev;
    if (estStDev <= nearZero) estStDev = nearZero;
    pzRow->calcMedian = calcMedian;
    pzRow->calcGSD = calcGSD;
    pzRow->scaledResidual = sqrt(size)*(pData[i].response - estMean) / estStDev;
    //double temp = bConstVar ? sqrt(sp / size) : pData[i].col4 / sqrt(size);
    double temp = pData[i].col4 / sqrt(size);
    pzRow->ebLower = calcMedian + gsl_cdf_tdist_Pinv(0.025, size - 1) * temp;
    pzRow->ebUpper = calcMedian + gsl_cdf_tdist_Pinv(0.975, size - 1) * temp;

#if defined(_DEBUG) || defined(DEBUGLOG)
    file << i << "  " << pData[i].dose << '\t' << size << '\t'
      << pData[i].response << '\t' << pData[i].col4 << '\t'
      << estMean << '\t' << estStDev << '\t' << pzRow->scaledResidual
      << '\t' << pzRow->ebLower << '\t' << pzRow->ebUpper << endl;
#endif

  } // for

#if defined(_DEBUG) || defined(DEBUGLOG)
  file << endl << endl;
  file.close();
#endif

  return nGroups;
}

// Returns number of estimated parameters for the model.
// - It requires the untransformed MLE values
int cModel::modelDF(Eigen::MatrixXd mleRaw, std::vector<bool> bFixed, bool *zBounded)
{
  //cout << "modelDF: initial nParms= " << nParms << endl;
  int estParmCount = nParms;
  for (int i = 0; i < nParms; i++) {
    zBounded[i] = false;
    double v = mleRaw(i, 0);
    if (bFixed[i] || bUnused[i]
        || fabs(v - zpriors(i, 3)) < BMDS_EPS
        || fabs(zpriors(i, 4) - v) < BMDS_EPS) {
      //cout << "modelDF: parm# " << i << " is fixed or bounded" << endl;
      --estParmCount;
      zBounded[i] = true;
    } // end if
  } // end for
  return estParmCount;
} // modelDF

// *********************  EXP2  MODEL    ****************************
bmds_exp2::bmds_exp2(bool bCV, bool bLN)
{
  nParms = nParmsBase;
  //cout << "model object base nparms= " << nParmsBase << endl;
  bConstVar = bCV;
  bLognormal = bLN;
  if (!bCV) {
    nParms++; // include rho
    //cout << "\tadding rho parm becuae not CV" << endl;
  }
  bUnused = new bool[nParms];
  // c and d are not used
  for (int i = 0; i < nParms; i++) bUnused[i] = false;
  bUnused[2] = bUnused[3] = true;
  //cout << "model object nparms= " << nParms << endl;
  return;
} // bmds_exp2(bool bCV, bool bLN)

bmds_exp2::~bmds_exp2() {}

void bmds_exp2::setParms(double * parms, bool bAdverseUp)
{
  a = parms[0];
  b = parms[1];
  c = parms[2];
  d = parms[3];
  sign = bAdverseUp ? 1 : -1;
  lnalpha = parms[nParms - 1];
  if (bConstVar) rho = 0;
  else rho = parms[nParms - 2];
} // bmds_exp2::setParms

// Model 2:     Y[dose] = a * exp{ sign * b * dose }
double bmds_exp2::mean(double dose)
{
  double result = a;
  if (dose > 0.0) result = a * exp(sign * b * dose);
  return result;
} // bmds_exp2::mean

// *********************  EXP3  MODEL    ****************************
bmds_exp3::bmds_exp3(bool bCV, bool bLN)
{
  nParms = nParmsBase;
  //cout << "model object base nparms= " << nParmsBase << endl;
  bConstVar = bCV;
  bLognormal = bLN;
  if (!bCV) {
    nParms++; // include rho
    //cout << "\tadding rho parm becuae not CV" << endl;
  }
  bUnused = new bool[nParms];
  // c is not used
  for (int i = 0; i < nParms; i++) bUnused[i] = false;
  bUnused[2] = true;
  //cout << "model object nparms= " << nParms << endl;
  return;
} // bmds_exp3(bool bCV, bool bLN)

bmds_exp3::~bmds_exp3() {}

void bmds_exp3::setParms(double * parms, bool bAdverseUp)
{
  a = parms[0];
  b = parms[1];
  c = parms[2];
  d = parms[3];
  sign = bAdverseUp ? 1 : -1;
  lnalpha = parms[nParms - 1];
  if (bConstVar) rho = 0;
  else rho = parms[nParms - 2];
} // bmds_exp3::setParms

// Model 3 : Y[dose] = a * exp{ sign * (b * dose) ^ d }
double bmds_exp3::mean(double dose)
{
  double result = a;
  if (dose > 0.0) result = a * exp(sign * pow((b * dose), d));
  return result;
} // bmds_exp3::mean

// *********************  EXP4  MODEL    ****************************
bmds_exp4::bmds_exp4(bool bCV, bool bLN)
{
  nParms = nParmsBase;
  //cout << "model object base nparms= " << nParmsBase << endl;
  bConstVar = bCV;
  bLognormal = bLN;
  if (!bCV) {
    nParms++; // include rho
    //cout << "\tadding rho parm becuae not CV" << endl;
  }
  bUnused = new bool[nParms];
  // d is not used
  for (int i = 0; i < nParms; i++) bUnused[i] = false;
  bUnused[3] = true;
  //cout << "model object nparms= " << nParms << endl;
  return;
} // bmds_exp4(bool bCV, bool bLN)

bmds_exp4::~bmds_exp4() {}

void bmds_exp4::setParms(double * parms, bool bAdverseUp)
{
  a = parms[0];
  b = parms[1];
  c = parms[2];
  d = parms[3];
  sign = bAdverseUp ? 1 : -1;
  lnalpha = parms[nParms - 1];
  if (bConstVar) rho = 0;
  else rho = parms[nParms - 2];
} // bmds_exp4::setParms

// Model 4 : Y[dose] = a * [c - (c - 1) * exp { -b * dose }]
double bmds_exp4::mean(double dose)
{
  double result = a;
  if (dose > 0.0) result = a * (c - (c - 1) * exp(-b * dose));
  return result;
} // bmds_exp4::mean

// *********************  EXP5  MODEL    ****************************
bmds_exp5::bmds_exp5(bool bCV, bool bLN)
{
  nParms = nParmsBase;
  //cout << "model object base nparms= " << nParmsBase << endl;
  bConstVar = bCV;
  bLognormal = bLN;
  if (!bCV) {
    nParms++; // include rho
    //cout << "\tadding rho parm becuae not CV" << endl;
  }
  bUnused = new bool[nParms];
  // all parameters are used
  for (int i = 0; i < nParms; i++) bUnused[i] = false;
  //cout << "model object nparms= " << nParms << endl;
  return;
} // bmds_exp5(bool bCV, bool bLN)

bmds_exp5::~bmds_exp5() {}

void bmds_exp5::setParms(double * parms, bool bAdverseUp)
{
  a = parms[0];
  b = parms[1];
  c = parms[2];
  d = parms[3];
  sign = bAdverseUp ? 1 : -1;
  lnalpha = parms[nParms - 1];
  if (bConstVar) rho = 0;
  else rho = parms[nParms - 2];
} // bmds_exp5::setParms

// Model 5 : Y[dose] = a * [c - (c - 1) * exp { -(b * dose) ^ d }]
double bmds_exp5::mean(double dose)
{
  cout << "exp5_mean: dose= " << dose << " a= " << a << " a= " << a << " b= "
    << b << " c= " << c << " d= " << d << endl;
  double result = a;
  if (dose > 0.0) result = a * (c - (c - 1) * exp(-pow((b*dose), d)));
  return result;
} // bmds_exp5::mean

// *********************  HILL  MODEL    ****************************
bmds_hill::bmds_hill(bool bCV, bool bLN)
{
  nParms = nParmsBase;
  //cout << "model object bConstVar= " << bCV << endl;
  bConstVar = bCV;
  bLognormal = bLN;
  if (!bCV) {
    nParms++; // include rho
    //cout << "\tadding rho parm becuae not CV" << endl;
  }
  bUnused = new bool[nParms];
  // all parameters are used
  for (int i = 0; i < nParms; i++) bUnused[i] = false;
  //cout << "model object nparms= " << nParms << endl;
  return;
} // bmds_hill(bool bCV, bool bLN)

bmds_hill::~bmds_hill() {}

void bmds_hill::setParms(double * parms, bool bAdverseUp)
{
  g = parms[0];
  v = parms[1];
  k = parms[2];
  n = parms[3];

  // lnalpha is always passed in via SetParms, regardless of model or variance model
  lnalpha = parms[nParms - 1];
  if (bConstVar) rho = 0;
  else rho = parms[nParms - 2];
} // bmds_hill::setParms

double bmds_hill::mean(double dose)
{
  double result = g;
  if (dose > 0.0) result = g + v * pow(dose, n) / (pow(k, n) + pow(dose, n));
  return result;
} // bmds_hill::mean

// *********************  Polynomial  MODEL    ****************************
bmds_poly::bmds_poly(int val, bool bCV, bool bLN)
{
  degree = val;
  nParms = nParmsBase + val + 1;
  //cout << "model object base nparms= " << nParmsBase << endl;
  bConstVar = bCV;
  bLognormal = bLN;
  if (!bCV) {
    nParms++; // include rho
    //cout << "\tadding rho parm becuae not CV" << endl;
  }
  p = new double[nParms];
  bUnused = new bool[nParms];
  // all parameters are used
  for (int i = 0; i < nParms; i++) bUnused[i] = false;
  //cout << "model object nparms= " << nParms << endl;
  return;
} // bmds_poly(int val, bool bCV, bool bLN)

bmds_poly::~bmds_poly()
{
  delete[] p;
}

void bmds_poly::setParms(double * parms, bool bAdverseUp)
{
  for (int i = 0; i < nParms; i++) {
    p[i] = parms[i];
  }

  // lnalpha is always passed in via SetParms, regardless of model or variance model
  lnalpha = parms[nParms - 1];
  if (bConstVar) rho = 0;
  else rho = parms[nParms - 2];
} // bmds_poly::setParms

// Y[dose] = beta_0 + beta_1*dose + beta_2*dose^2 + ...
double bmds_poly::mean(double dose)
{
  double result = p[degree];
  for (int i = degree - 1; i >= 0; i--) {
    result = result * dose + p[i];
  } // end for
  return result;
} // bmds_poly::mean

// *********************  Power  MODEL    ****************************
bmds_power::bmds_power(bool bCV, bool bLN)
{
  nParms = nParmsBase;
  //cout << "model object base nparms= " << nParmsBase << endl;
  bConstVar = bCV;
  bLognormal = bLN;
  if (!bCV) {
    nParms++; // include rho
    //cout << "\tadding rho parm becuae not CV" << endl;
  }
  bUnused = new bool[nParms];
  // all parameters are used
  for (int i = 0; i < nParms; i++) bUnused[i] = false;
  //cout << "model object nparms= " << nParms << endl;
  return;
} // bmds_power(bool bCV, bool bLN)

bmds_power::~bmds_power() {}

void bmds_power::setParms(double * parms, bool bAdverseUp)
{
  g = parms[0];
  b = parms[1];
  n = parms[2];

  // lnalpha is always passed in via SetParms, regardless of model or variance model
  lnalpha = parms[nParms - 1];
  if (bConstVar) rho = 0;
  else rho = parms[nParms - 2];
} // bmds_power::setParms

// Y[dose] = control + slope * dose^power
double bmds_power::mean(double dose)
{
  double result = g;
  if (dose > 0.0) result = g + b * pow(dose, n);
  return result;
} // bmds_power::mean

