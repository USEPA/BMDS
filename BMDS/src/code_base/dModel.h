#pragma once

#ifndef DMODELH
#define DMODELH


#include "bmds_entry.h"
#include "bmds_dmodels.h"


class dModel {
public:
  dModel();
  ~dModel();
  virtual void setParms(double *parms) =0;
  virtual double mean(double dose) =0;
  void calcGoF(dGoF_t *zOut, BMDSInputData_t *dataIn, const int nRows, int npEst);

private:
  const int nParms = 0;
};

class bmds_loglogistic: public dModel {
public:
  bmds_loglogistic();
  ~bmds_loglogistic();
  void setParms(double *parms);
  double mean(double dose);

private:
  const int nParms = 3;
  double g, a, b;
};

class bmds_gamma : public dModel {
public:
  bmds_gamma();
  ~bmds_gamma();
  void setParms(double *parms);
  double mean(double dose);

private:
  const int nParms = 3;
  double g, a, b;
};

class bmds_dhill : public dModel {
public:
  bmds_dhill();
  ~bmds_dhill();
  void setParms(double *parms);
  double mean(double dose);

private:
  const int nParms = 4;
  double g, n, a, b;
};

class bmds_logistic : public dModel {
public:
  bmds_logistic();
  ~bmds_logistic();
  void setParms(double *parms);
  double mean(double dose);

private:
  const int nParms = 2;
  double a, b;
};

class bmds_logprobit : public dModel {
public:
  bmds_logprobit();
  ~bmds_logprobit();
  void setParms(double *parms);
  double mean(double dose);

private:
  const int nParms = 3;
  double g, a, b;
};

class bmds_probit : public dModel {
public:
  bmds_probit();
  ~bmds_probit();
  void setParms(double *parms);
  double mean(double dose);

private:
  const int nParms = 2;
  double a, b;
};

class bmds_qlinear : public dModel {
public:
  bmds_qlinear();
  ~bmds_qlinear();
  void setParms(double *parms);
  double mean(double dose);

private:
  const int nParms = 2;
  double g, b;
};

class bmds_weibull : public dModel {
public:
  bmds_weibull();
  ~bmds_weibull();
  void setParms(double *parms);
  double mean(double dose);

private:
  const int nParms = 3;
  double g, a, b;
};

class bmds_multistage : public dModel {
public:
  bmds_multistage(int degree);
  ~bmds_multistage();
  void setParms(double *parms);
  double mean(double dose);

private:
  int nParms = 0; // Depends on degree
  int degree = 0; // Set by the constructor
  double *p;
};

#endif
