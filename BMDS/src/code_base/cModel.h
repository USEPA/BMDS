#pragma once
#include <vector>

class cModel {
public:
  cModel();
  ~cModel();
  virtual void setParms(double *parms, bool bAdverseUp) = 0;
  virtual double mean(double dose) = 0;
  void setPriors(Eigen::MatrixXd priors);
  int calcGoF(cGoFRow_t *zOut, BMDSInputType_t inputType, BMDSInputData_t *dataIn,
               const int nRows, int npEst, Eigen::MatrixXd X, Eigen::MatrixXd Y);
  int modelDF(Eigen::MatrixXd mleRaw, std::vector<bool> bFixed, bool *zBounded);

protected:
  // To maintain BMDS 2.x compatibility, alpha is used for constant variance, and
  // lnalpha is used for modeled variance. EXP models always use lnalpha regardless.
  // -1 is not a valid value for alpha; so we use it as a dummy flag value indicating
  // that lnalpha should be used.
  double rho, alpha, lnalpha;
  int nParms = 0;
  bool bConstVar;
  bool bLognormal;
  bool *bUnused; // Indicates whether model parms are used
  Eigen::MatrixXd zpriors;

private:
};

class bmds_exp2 : public cModel {
public:
  bmds_exp2(bool bConstVar, bool bLognormal);
  ~bmds_exp2();
  void setParms(double *parms, bool bAdverseUp);
  double mean(double dose);

private:
  double a, b, c, d;
  double sign; // based on direction of adversity
  // base count includes alpha parameter
  const int nParmsBase = 5;
};

class bmds_exp3 : public cModel {
public:
  bmds_exp3(bool bConstVar, bool bLognormal);
  ~bmds_exp3();
  void setParms(double *parms, bool bAdverseUp);
  double mean(double dose);

private:
  double a, b, c, d;
  double sign; // based on direction of adversity
  // base count includes alpha parameter
  const int nParmsBase = 5;
};

class bmds_exp4 : public cModel {
public:
  bmds_exp4(bool bConstVar, bool bLognormal);
  ~bmds_exp4();
  void setParms(double *parms, bool bAdverseUp);
  double mean(double dose);

private:
  double a, b, c, d;
  double sign; // based on direction of adversity
// base count includes alpha parameter
  const int nParmsBase = 5;
};

class bmds_exp5 : public cModel {
public:
  bmds_exp5(bool bConstVar, bool bLognormal);
  ~bmds_exp5();
  void setParms(double *parms, bool bAdverseUp);
  double mean(double dose);

private:
  double a, b, c, d;
  double sign; // based on direction of adversity
// base count includes alpha parameter
  const int nParmsBase = 5;
};

class bmds_hill : public cModel {
public:
  bmds_hill(bool bConstVar, bool bLognormal);
  ~bmds_hill();
  void setParms(double *parms, bool bAdverseUp);
  double mean(double dose);

private:
  double g, v, k, n;
  // base count includes alpha parameter
  const int nParmsBase = 5;
};

class bmds_poly : public cModel {
public:
  bmds_poly(int degree, bool bConstVar, bool bLognormal);
  ~bmds_poly();
  void setParms(double *parms, bool bAdverseUp);
  double mean(double dose);

private:
  // base count includes alpha parameter
  const int nParmsBase = 1; // Depends on degree
  int degree = 0; // Set by the constructor
  double *p; // parameter value array
};

class bmds_power : public cModel {
public:
  bmds_power(bool bConstVar, bool bLognormal);
  ~bmds_power();
  void setParms(double *parms, bool bAdverseUp);
  double mean(double dose);

private:
  double g, b, n;
  // base count includes alpha parameter
  const int nParmsBase = 4;
};

