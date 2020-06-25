#pragma once
#ifndef bmds_dmodelsH
#define bmds_dmodelsH

#include "bmd_calculate.h"
#ifdef R_COMPILATION
#define LOGFILENAME "BOBSUCKS.OUT" 
#else
#define LOGFILENAME __FUNCTION__ ## ".log"
#endif // R_COMPILATION

#define NPARM_DHILL 4
#define NPARM_GAMMA 3
#define NPARM_LNLOGISTIC 3
#define NPARM_LNPROBIT 3
#define NPARM_LOGISTIC 2
#define NPARM_MSTAGE_BASE 1 // Base # of parms exluding betas according to degree
#define NPARM_MSTAGE(x) (1 + (x))
#define NPARM_PROBIT 2
#define NPARM_QLINEAR 2
#define NPARM_WEIBULL 3

#define iBACKGROUND 0 // Array index for background parameter for all models
//constexpr auto BMDS_EPS = 1.523e-8;
constexpr auto BMDS_EPS = 1.0e-6;
constexpr auto BMDS_ONE = 1 - BMDS_EPS;

static const char* MODELNAME[] = { NULL, "D-hill", "Gamma", "Logistic", "Log-logistic",
                                  "Log-probit", "Multistage", "Probit", "Quantal-linear",
                                  "Weibull" };

int modelDF(bmd_analysis *fitout, Eigen::MatrixXd priors, vector<bool>& bFixed, bool *zBounded);
double fixedParmToLogistic(double x);
void analysisOfDeviance(DichotomousDeviance_t *zOut, Eigen::MatrixXd& Y,
                        Eigen::MatrixXd& X, double llFit, int nEstParms);

#endif // bmds_dmodelsH
