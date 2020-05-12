#pragma once
#include "bmds_entry.h"

class cInputData {
public:
  cInputData(BMDSInputData_t *dataIn, BMDSInputType_t inputType, int n, bool bLN);
  ~cInputData();

//protected:
  double doseMin = DBL_MAX;
  double doseMax = 0;
  double yScale = 1;  // Value used to scale the responses
  int nDoses = 0;    // Number of dose groups
  BMDSInputType_t origType;
  bool bLognormal; // Lognormal distribution is assumed
  bool bScaleResponse = true; // Flag that controls whether responses are scaled
  Eigen::MatrixXd origData;
  Eigen::MatrixXd drData;

private:
  void importSummaryData(BMDSInputData_t *dataIn, int n, bool bLN);
  void importIndividualData(BMDSInputData_t *dataIn, int n, bool bLN);
};