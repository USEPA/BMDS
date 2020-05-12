#include "stdafx.h"
#ifdef R_COMPILATION
#  include <RcppEigen.h>
#  include <RcppGSL.h>
#else 
#  include <Eigen/Dense>
#endif // R_COMPILATION
#include "gsl/gsl_errno.h"
#include "gsl/gsl_cdf.h"
#include <vector>
#include <iostream>
#include "cInputData.h"

using namespace std;
using namespace Eigen;

#define LOGFILENAME = __FUNCTION__ ## ".log"

// Indices for imported d-r data
constexpr auto DOSE = 0; // Summarized and individual responses
constexpr auto MEAN = 1; // Summarized responses
constexpr auto RESP = 1; // Individual responses
constexpr auto STDV = 2; // Summarized responses
constexpr auto LNRS = 2; // Individual responses
constexpr auto SZCO = 3; // Summarized responses

void createSummaryData(Eigen::MatrixXd& orig, Eigen::MatrixXd& summary);

cInputData::cInputData(BMDSInputData_t *dataIn, BMDSInputType_t inputType, int n, bool bLN)
{
  //int nDoses = 0;    // Number of dose groups
  //double doseMin;
  //double doseMax;
  //double yScale;  // Value used to scale the responses

  bLognormal = bLN;
  origType = inputType;
  bool bSS = (inputType == eCont_4);

  // Currently, assume that dose groups are sorted, but we might
  // want to sort explicitly to be safe!!
  cout << "bSS= " << bSS << endl << flush;
  if (bSS) importSummaryData(dataIn, n, bLN);
  else importIndividualData(dataIn, n, bLN);
  nDoses = drData.rows();
} // cInputData()

cInputData::~cInputData() {}

// Copy the input data into a matrix and get the min/max doses at the same time.
void cInputData::importSummaryData(BMDSInputData_t *dataIn, int n, bool bLN)
{
  origData.resize(n, 4);
  drData.resize(n, 4);
  for (int i = 0; i < n; i++) {
    double dose = dataIn[i].dose;
    if (dose > doseMax) doseMax = dose;
    else if (dose < doseMin) doseMin = dose;
    origData(i, DOSE) = dose;
    origData(i, MEAN) = dataIn[i].response; // Group mean
    origData(i, STDV) = dataIn[i].col4;     // Std dev
    origData(i, SZCO) = dataIn[i].groupSize;
  } // for i
  if (bScaleResponse) yScale = origData(0, MEAN); // Scale responses by mean for the lowest dose
  drData.col(DOSE) = origData.col(DOSE) / doseMax;
  drData.col(MEAN) = origData.col(MEAN) / yScale;
  drData.col(STDV) = origData.col(STDV) / yScale;
  drData.col(SZCO) = origData.col(SZCO);
  // If lognormal distribution, approximate responses on the log scale
  if (bLN) {
    for (int i = 0; i < n; i++) {
      double x = drData(i, MEAN);
      double y = drData(i, STDV);
      drData(i, MEAN) = log(x) - log(1 + pow(y / x, 2.0)) / 2;
      drData(i, STDV) = sqrt(log(1.0 + pow(y / x, 2.0)));
    } // for i
  } // if (bLN)

  return;
} // importSummaryData()

// Copy individual d-r data into a matrix
// - Get the min and max doses
// - If lognormal distribution, convert to log scale
void cInputData::importIndividualData(BMDSInputData_t *dataIn, int n, bool bLN)
{
  origData.resize(n, 2);
  for (int i = 0; i < n; i++) {
    double dose = dataIn[i].dose;
    if (dose > doseMax) doseMax = dose;
    else if (dose < doseMin) doseMin = dose;
    origData(i, DOSE) = dose;
    origData(i, RESP) = dataIn[i].response;
  } // for i
  if (bLN) origData.col(RESP).array() = origData.col(RESP).array().log();
  createSummaryData(origData, drData);
  if (bScaleResponse) yScale = drData(0, MEAN); // Scale responses by mean for the lowest dose

  drData.col(DOSE) = drData.col(DOSE) / doseMax;
  drData.col(MEAN) = drData.col(MEAN) / yScale;
  drData.col(STDV) = drData.col(STDV) / yScale;


  return;
}

void createSummaryData(Eigen::MatrixXd& orig, Eigen::MatrixXd& summary)
{
  int dg = 0; // Rows in the summarized data set
  // Get the number of dose groups
  // Assume that doses are sorted
  double dose = orig(0, DOSE);
  dg = 1;
  for (int i = 1; i < orig.rows(); i++) {
    if (orig(i, DOSE) != dose) {
      // New dose group
      dose = orig(i, DOSE);
      ++dg;
    } // if
  } // for i

  summary.resize(dg, 4);
  dose = orig(0, DOSE);
  int gSize = 1; // Number of observations in a dose group
  double sum1 = orig(0, RESP); // sum of responses in group
  double sum2 = orig(0, RESP) * orig(0, RESP); // sum squares of responses
  dg = 0;
  // NOTE: Loop starts at the 2nd element!!!!
  for (int i = 1; i < orig.rows(); i++) {
    if (orig(i, DOSE) != dose) { // New dose group
      // Wrap up the current group
      summary(dg, DOSE) = dose;
      summary(dg, MEAN) = sum1 / gSize;
      summary(dg, SZCO) = gSize;
      summary(dg, STDV) = sqrt((sum2 - sum1 * sum1 / gSize) / (gSize - 1));
      // Start a new dose group
      ++dg;
      dose = orig(i, DOSE);
      gSize = 0;
      sum1 = sum2 = 0;
    } // if
    ++gSize;
    sum1 += orig(i, RESP);
    sum2 += orig(i, RESP) * orig(i, RESP);
  } // for
  // Wrap up the last dose group
  summary(dg, DOSE) = dose;
  summary(dg, MEAN) = sum1 / gSize;
  summary(dg, SZCO) = gSize;
  summary(dg, STDV) = sqrt((sum2 - sum1 * sum1 / gSize) / (gSize - 1));


  return;
} // createSummaryData()
