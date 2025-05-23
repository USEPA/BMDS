#ifdef WIN32
#  include "pch.h"
#else
#  include "stdafx.h"
#endif  // Precompiled header - does nothing if building R version
#include "log_likelihoods.h"
#include "lognormalModels.h"
#include "lognormal_POLYNOMIAL_NC.h"

#ifdef R_COMPILATION
// necessary things to run in R
#  include <RcppEigen.h>
#  include <RcppGSL.h>
#else
#  include <Eigen/Dense>
#endif

#include <gsl/gsl_randist.h>

#include "gsl/gsl_cdf.h"

/////////////////////////////////////////////////////////////////////////
// function: type_of_profile()
// purpose: Tell the optimizer which profile likelihood method is
//			best for the given bmdtype.  For all models, the HYBRID
//          is always represented as an equality constraint.
// input:
//	contbmd TYPE
// output:
//  PROFILE_INEQUALITY - One of the parameters can be made equal to the others as a function
//						 of the fixed BMD. The optimizer thus optimizes a
// smaller problem
//  PROFILE_EQUALITY   - The BMD is  a function of multiple parameters and can not be disentangled
//                       An equality constraint is used here.
//////////////////////////////////////////////////////////////////////////
int lognormalPOLYNOMIAL_BMD_NC::type_of_profile(contbmd TYPE) {
  switch (TYPE) {
    case CONTINUOUS_BMD_ABSOLUTE:
    case CONTINUOUS_BMD_STD_DEV:
    case CONTINUOUS_BMD_REL_DEV:
    case CONTINUOUS_BMD_POINT:
    case CONTINUOUS_BMD_EXTRA:
    case CONTINUOUS_BMD_HYBRID_EXTRA:
    default:
      return PROFILE_EQUALITY;  // HYBRID IS ALWAYS AN EQUALITY CONSTRAINT
      break;
  }
}

/////////////////////////////////////////////////////////////////////////
// function: parameter_to_remove()
// purpose: Tell the optimizer which profile likelihood method is
//			best for the given bmdtype.  For all models, the HYBRID
//          is always represented as an equality constraint.
// input:
//	contbmd TYPE
// output:
//  PROFILE_INEQUALITY - One of the parameters can be made equal to the others as a function
//						 of the fixed BMD. The optimizer thus optimizes a
// smaller problem
//  PROFILE_EQUALITY   - The BMD is  a function of multiple parameters and can not be disentangled
//                       An equality constraint is used here.
//////////////////////////////////////////////////////////////////////////
int lognormalPOLYNOMIAL_BMD_NC::parameter_to_remove(contbmd TYPE) {
  switch (TYPE) {
    case CONTINUOUS_BMD_ABSOLUTE:
    case CONTINUOUS_BMD_REL_DEV:
      return 1;  // beta parameter
    case CONTINUOUS_BMD_POINT:
    case CONTINUOUS_BMD_EXTRA:
      return 0;
      break;
    case CONTINUOUS_BMD_STD_DEV:
      return nParms() - 1;
      break;
    case CONTINUOUS_BMD_HYBRID_EXTRA:
    default:
      return -1;  // NO PARAMETER IS REMOVED THUS IT IS A NEGATIVE
                  // INDEX AND WILL THROW AN ERROR
      break;
  }
}

/*Function: bmd_start_extra_absolute
* Purpose : Give a starting value that satisfied the equality constraints so a profile likelihood
*           can be performed. This function is used in an optimization that finds the
                        closest point to the supplied value (usually the previous MLE)
                        that satisfies the equality constraint. Note it always modifies parameters
that are essentially unbounded - nu
*/
double lognormalPOLYNOMIAL_BMD_NC::bmd_start_absolute(
    unsigned n, const double *b, double *grad, void *data
) {
  start_data *sdata = (start_data *)data;
  Eigen::MatrixXd theta = sdata->theta;

  double returnV = pow(theta(deg + 1, 0) - b[deg + 1], 2.0) + pow(theta(0, 0) - b[0], 2.0);
  if (!sdata->isIncreasing) sdata->BMRF *= -1;

  double temp = 0.0;
  double sum = 0.0;

  if (deg == 1) {
    temp = sdata->BMRF / sdata->BMD;
  } else {
    for (int i = 2; i < deg + 1; i++) {
      sum += b[i] * pow(sdata->BMD, i);
      returnV += pow(theta(i, 0) - b[i], 2.0);
    }

    temp = (sdata->BMRF - sum) / sdata->BMD;
  }

  returnV += pow(temp - theta(1, 0), 2.0);

  return returnV;
}

std::vector<double> lognormalPOLYNOMIAL_BMD_NC::bmd_start_absolute_clean(
    std::vector<double> x, double BMRF, double BMD, bool isIncreasing
) {
  if (!isIncreasing) BMRF *= -1.0;

  double temp = 0.0;
  double sum = 0.0;

  if (deg == 1) {
    temp = BMRF / BMD;
  } else {
    for (int i = 2; i < deg + 1; i++) {
      sum += x[i] * pow(BMD, i);
    }

    temp = (BMRF - sum) / BMD;
  }
  x[1] = temp;
  return x;
}

/******************************************
* Function: bmd_start_reldev
* Purpose : Give a starting value that satisfied the equality constraints so a profile likelihood
*           can be performed. This function is used in an optimization that finds the
closest point to the supplied value (usually the previous MLE)
that satisfies the equality constraint. Note it always modifies parameters that are essentially
unbounded - nu
*/
double lognormalPOLYNOMIAL_BMD_NC::bmd_start_reldev(
    unsigned n, const double *b, double *grad, void *data
) {
  /*	double t;
          start_data 		 *sdata = (start_data *)data;

          if (sdata->isIncreasing)
                  t = sdata->BMRF;
          else
                  t = 1.0 - sdata->BMRF;
          // key off of mu in this one

          Eigen::MatrixXd  theta = sdata->theta;


          Eigen::MatrixXd theta_2 = theta;
          for (int i = 0; i < n; i++) { theta_2(i, 0) = b[i]; }

          double temp = 0.0;
          for (int i = 2; i < theta.rows(); i++)
                  temp += b[i] * pow(sdata->BMD, i);

          temp = sdata->isIncreasing ? (t * b[0] - temp) : (-t * b[0] - temp);
          temp = temp / (sdata->BMD);

          theta_2(1, 0) = temp;
          double returnV = 0.0;
          for (int i = 0; i < theta.rows(); i++)
                  returnV += pow(theta(i, 0) - theta_2(i, 0), 2.0);


          return  returnV;
  */
  double t;
  start_data *sdata = (start_data *)data;

  if (sdata->isIncreasing)
    t = sdata->BMRF;
  else
    t = 1.0 - sdata->BMRF;
  // key off of mu in this one

  Eigen::MatrixXd theta = sdata->theta;

  Eigen::MatrixXd theta_2 = theta;
  for (int i = 0; i < n; i++) {
    theta_2(i, 0) = b[i];
  }

  double temp = 0.0;

  // int deg_poly = sdata->M->isConstVar():theta.rows():theta.rows()-1;

  for (int i = 1; i <= deg; i++) temp += b[i] * pow(sdata->BMD, i);

  temp = sdata->isIncreasing ? temp / t : -temp / t;

  theta_2(0, 0) = temp;
  double returnV = 0.0;

  // Euclidean distance between the values
  for (int i = 0; i < theta.rows(); i++) returnV += pow(theta(i, 0) - theta_2(i, 0), 2.0);

  // cout << theta << endl;
  return returnV;
}

std::vector<double> lognormalPOLYNOMIAL_BMD_NC::bmd_start_reldev_clean(
    std::vector<double> x, double BMRF, double BMD, bool isIncreasing
) {
  /*double t;

  Eigen::MatrixXd theta_2(x.size(), 1);
  for (int i = 0; i < x.size(); i++) { theta_2(i, 0) = x[i]; }
  Eigen::MatrixXd d(2, 1); d << 0.0, BMD;
  Eigen::MatrixXd mu = mean(theta_2, d);

  if (isIncreasing)
          t = BMRF;
  else
          t = 1.0 - BMRF;

  double temp = 0.0;
  for (int i = 2; i < theta_2.rows(); i++)
          temp += x[i] * pow(BMD, i);
  temp = isIncreasing ? (t * x[0] - temp) : (-temp - t * x[0]);
  temp = temp / (BMD);


  x[1] = temp;

  return x;*/
  double t;

  Eigen::MatrixXd theta_2(x.size(), 1);
  for (int i = 0; i < x.size(); i++) {
    theta_2(i, 0) = x[i];
  }
  Eigen::MatrixXd d(2, 1);
  d << 0.0, BMD;
  Eigen::MatrixXd mu = mean(theta_2, d);

  if (isIncreasing)
    t = BMRF;
  else
    t = 1.0 - BMRF;

  double temp = 0.0;
  // loops over the degree of the polynomial
  for (int i = 1; i <= deg; i++) temp += x[i] * pow(BMD, i);

  temp = isIncreasing ? (temp / t) : -(temp / t);

  x[0] = temp;

  return x;
}

/************************************************************
* Function: bmd_start_stddev
* Purpose : Give a starting value that satisfied the equality constraints so a profile likelihood
*           can be performed. This function is used in an optimization that finds the
closest point to the supplied value (usually the previous MLE)
that satisfies the equality constraint. Note it always modifies parameters that are essentially
unbounded - nu
**************************************************************/
double lognormalPOLYNOMIAL_BMD_NC::bmd_start_stddev(
    unsigned n, const double *b, double *grad, void *data
) {
  // key off of mu in this one
  start_data *sdata = (start_data *)data;
  Eigen::MatrixXd theta = sdata->theta;
  if (!sdata->isIncreasing) sdata->BMRF *= -1.0;

  Eigen::MatrixXd d(2, 1);
  d << 0.0, sdata->BMD;
  Eigen::MatrixXd theta_2 = theta;
  for (int i = 0; i < n; i++) {
    theta_2(i, 0) = b[i];
  }
  Eigen::MatrixXd mu = mean(theta_2, d);
  mu = exp(mu.array());

  //	Eigen::MatrixXd md = fabs(exp(log(mu(0,0)) + BMRF * pow(t, 0.5)) - mu(0,0));
  double temp = fabs(mu(1, 0) - mu(0, 0)) / mu(0, 0);
  if (sdata->isIncreasing)
    temp = temp + 1;
  else
    temp = 1 - temp;
  ;

  temp = 2.0 * log(log(temp) / sdata->BMRF);

  double returnV = 0.0;

  // compute the squared Euclidean Distance
  returnV += pow(temp - theta(n - 1, 0), 2.0);  // variance parameter
  for (int i = 0; i < n - 1; i++) {
    returnV += pow(theta(i, 0) - b[i], 2.0);
  }

  return returnV;
}

std::vector<double> lognormalPOLYNOMIAL_BMD_NC::bmd_start_stddev_clean(
    std::vector<double> x, double BMRF, double BMD, bool isIncreasing
) {
  if (!isIncreasing) BMRF *= -1.0;

  // key off of the std dev in this one
  Eigen::MatrixXd d(2, 1);
  d << 0.0, BMD;
  Eigen::MatrixXd theta_2(x.size(), 1);
  ;
  for (int i = 0; i < x.size(); i++) {
    theta_2(i, 0) = x[i];
  }
  Eigen::MatrixXd mu = mean(theta_2, d);
  mu = exp(mu.array());

  double temp = fabs(mu(1, 0) - mu(0, 0)) / mu(0, 0);
  if (isIncreasing)
    temp = temp + 1;
  else
    temp = 1 - temp;

  temp = 2.0 * log(log(temp) / BMRF);

  x[x.size() - 1] = temp;
  return x;
}

/***********************************************
* Function: bmd_start_extra_relative
* Purpose : Give a starting value that satisfied the equality constraints so a profile likelihood
*           can be performed. This function is used in an optimization that finds the
closest point to the supplied value (usually the previous MLE)
that satisfies the equality constraint. Note it always modifies parameters that are essentially
unbounded - nu
*/
double lognormalPOLYNOMIAL_BMD_NC::bmd_start_point(
    unsigned n, const double *b, double *grad, void *data
) {
  // key off background in this one
  start_data *sdata = (start_data *)data;
  Eigen::MatrixXd theta = sdata->theta;

  double returnV = pow(theta(deg + 1, 0) - b[deg + 1], 2.0);

  double temp = 0.0;
  double sum = 0.0;

  for (int i = 1; i < deg + 1; i++) {
    sum += b[i] * pow(sdata->BMD, i);
    returnV += pow(theta(i, 0) - b[i], 2.0);
  }

  temp = sdata->BMRF - sum;
  returnV += pow(temp - theta(0, 0), 2.0);

  return returnV;
}

std::vector<double> lognormalPOLYNOMIAL_BMD_NC::bmd_start_point_clean(
    std::vector<double> x, double BMRF, double BMD, bool isIncreasing
) {
  double temp = 0.0;
  double sum = 0.0;

  for (int i = 1; i < deg + 1; i++) {
    sum += x[i] * pow(BMD, i);
  }

  temp = BMRF - sum;
  x[0] = temp;
  return x;
}

/********************************************
* Function: bmd_start_extra
* Purpose : Give a starting value that satisfied the equality constraints so a profile likelihood
*           can be performed. This function is used in an optimization that finds the
closest point to the supplied value (usually the previous MLE)
that satisfies the equality constraint. Note it always modifies parameters that are essentially
unbounded - nu
*********************************************/
double lognormalPOLYNOMIAL_BMD_NC::bmd_start_extra(
    unsigned n, const double *b, double *grad, void *data
) {
  return 0.0;
}

std::vector<double> lognormalPOLYNOMIAL_BMD_NC::bmd_start_extra_clean(
    std::vector<double> x, double BMRF, double BMD, bool isIncreasing
) {
  return x;
}

/*******************************************************
*Function: bmd_start_hybrid_extra
* Purpose : Give a starting value that satisfied the equality constraints so a profile likelihood
*           can be performed. This function is used in an optimization that finds the
closest point to the supplied value (usually the previous MLE)
that satisfies the equality constraint. Note it always modifies parameters that are essentially
unbounded - nu
*********************************************************/
double lognormalPOLYNOMIAL_BMD_NC::bmd_start_hybrid_extra(
    unsigned n, const double *b, double *grad, void *data
) {
  start_data *sdata = (start_data *)data;
  double NOT_ADVERSE_P = 1.0 - sdata->tail_prob;
  double TAIL_PROB = sdata->tail_prob;
  Eigen::MatrixXd theta = sdata->theta;
  Eigen::MatrixXd theta2 = theta;
  /////////////////////////////////////////////////////////////////////////////////////
  for (int i = 0; i < n; i++) {
    theta2(i, 0) = b[i];
  }
  /////////////////////////////////////////////////////////////////////////////////////
  Eigen::MatrixXd d(2, 1);
  d << 0, sdata->BMD;
  Eigen::MatrixXd temp_mean = mean(theta2, d);
  Eigen::MatrixXd temp_var = variance(theta2, d);
  /////////////////////////////////////////////////////////////////////////////////////
  double mu_zero = temp_mean(0, 0);
  double std_zero = sqrt(temp_var(0, 0));
  double ct_off =
      gsl_cdf_lognormal_Pinv(sdata->isIncreasing ? NOT_ADVERSE_P : TAIL_PROB, mu_zero, std_zero);
  /////////////////////////////////////////////////////////////////////////////////////
  double returnV = 0.0;
  double temp;
  double k1 = gsl_cdf_ugaussian_Pinv(NOT_ADVERSE_P * sdata->BMRF + TAIL_PROB);
  double k0 = gsl_cdf_ugaussian_Pinv(TAIL_PROB);

  // Need to differenciate between increasing and decreasing
  if (sdata->isIncreasing) {
    temp = (temp_mean(1, 0) - temp_mean(0, 0)) / (k1 - k0);  // This is the standard deviation
    temp = 2.0 * log(temp);  // transform it to the log scale and make it a variance

  } else {
    temp = (temp_mean(1, 0) - temp_mean(0, 0)) / (k0 - k1);  // This is the standard deviation
    temp = 2.0 * log(temp);  // transform it to the log scale and make it a variance
  }

  //////////////////////////////////////////////////////////////////////
  for (int i = 0; i <= n - 2; i++) {
    returnV += pow(theta(i, 0) - b[i], 2.0);
  }
  returnV += pow(temp - theta(n - 1, 0), 2.0);

  return returnV;
}

std::vector<double> lognormalPOLYNOMIAL_BMD_NC::bmd_start_hybrid_extra_clean(
    std::vector<double> x, double BMRF, double BMD, bool isIncreasing, double tail_prob
) {
  double NOT_ADVERSE_P = 1.0 - tail_prob;
  double TAIL_PROB = tail_prob;
  Eigen::MatrixXd theta2(x.size(), 1);
  /////////////////////////////////////////////////////////////////////////////////////
  for (int i = 0; i < x.size(); i++) {
    theta2(i, 0) = x[i];
  }
  /////////////////////////////////////////////////////////////////////////////////////
  Eigen::MatrixXd d(2, 1);
  d << 0, BMD;
  Eigen::MatrixXd temp_mean = mean(theta2, d);
  Eigen::MatrixXd temp_var = variance(theta2, d);
  /////////////////////////////////////////////////////////////////////////////////////
  double mu_zero = temp_mean(0, 0);
  double std_zero = sqrt(temp_var(0, 0));
  double ct_off =
      gsl_cdf_lognormal_Pinv(isIncreasing ? NOT_ADVERSE_P : TAIL_PROB, mu_zero, std_zero);
  /////////////////////////////////////////////////////////////////////////////////////
  double returnV = 0.0;
  double temp;
  double k1 = gsl_cdf_ugaussian_Pinv(NOT_ADVERSE_P * BMRF + TAIL_PROB);

  double k0 = gsl_cdf_ugaussian_Pinv(TAIL_PROB);

  // Need to differenciate between increasing and decreasing
  if (isIncreasing) {
    temp = (temp_mean(1, 0) - temp_mean(0, 0)) / (k1 - k0);  // This is the standard deviation
    // cout << temp << endl;
    temp = 2.0 * log(temp);  // transform it to the log scale and make it a variance

  } else {
    temp = (temp_mean(1, 0) - temp_mean(0, 0)) / (k0 - k1);  // This is the standard deviation
    temp = 2.0 * log(temp);  // transform it to the log scale and make it a variance
  }

  //////////////////////////////////////////////////////////////////////
  x[x.size() - 1] = temp;  // last term is ALWAYS the variance
  return x;
}

/****************************************************************
// Functions: lognormalPOLYNOMIAL_BMD_NC::bmd_absolute(Eigen::MatrixXd theta, double BMRF, bool
isIncreasing)
//			  lognormalPOLYNOMIAL_BMD_NC::bmd_stdev(Eigen::MatrixXd theta, double BMRF,
bool isIncreasing)
//			  lognormalPOLYNOMIAL_BMD_NC::bmd_reldev(Eigen::MatrixXd theta, double BMRF,
bool isIncreasing)
//            lognormalPOLYNOMIAL_BMD_NC::bmd_point(Eigen::MatrixXd theta, double BMRF, bool
isIncreasing)
// Purpose :  return the BMD given the parameter values theta and the BMRF. Note they are  call the
code
//            in  lognormalPOLYNOMIAL_BMD_NC::bmd_absolute(Eigen::MatrixXd theta, double BMRF, bool
isIncreasing)
//
* ***************************************************************/
double lognormalPOLYNOMIAL_BMD_NC::bmd_absolute_bound(
    Eigen::MatrixXd theta, double BMD, double BMRF, bool isIncreasing
) {
  Eigen::MatrixXd d(2, 1);
  d << 0.0, BMD;
  Eigen::MatrixXd temp = mean(theta, d);
  temp = exp(temp.array());

  double rValue = fabs(temp(0, 0) - temp(1, 0)) - BMRF;
  return rValue;
}

double lognormalPOLYNOMIAL_BMD_NC::bmd_stdev_bound(
    Eigen::MatrixXd theta, double BMD, double BMRF, bool isIncreasing
) {
  if (!isIncreasing) {
    BMRF *= -1.0;
  }
  Eigen::MatrixXd d(1, 1);
  d(0, 0) = 0.0;
  Eigen::MatrixXd temp = variance(theta, d);
  Eigen::MatrixXd med = mean(theta, d);
  med = exp(med.array());
  double t = temp(0, 0);
  Eigen::MatrixXd md = abs(exp(log(med.array()) + BMRF * pow(t, 0.5)) - med.array());
  return bmd_absolute_bound(theta, BMD, md(0, 0), isIncreasing);
}

double lognormalPOLYNOMIAL_BMD_NC::bmd_reldev_bound(
    Eigen::MatrixXd theta, double BMD, double BMRF, bool isIncreasing
) {
  Eigen::MatrixXd d(1, 1);
  d(0, 0) = 0.0;
  Eigen::MatrixXd temp = mean(theta, d);
  temp = exp(temp.array());

  double t;

  if (isIncreasing)
    t = BMRF * temp(0, 0);
  else
    t = temp(0, 0) * (1.0 - BMRF);

  return bmd_absolute_bound(theta, BMD, t, isIncreasing);
}

double lognormalPOLYNOMIAL_BMD_NC::bmd_extra_bound(
    Eigen::MatrixXd theta, double BMD, double BMRF, bool isIncreasing
) {
  Eigen::MatrixXd d(1, 1);
  d(0, 0) = 0.0;
  Eigen::MatrixXd temp = mean(theta, d);
  temp = exp(temp.array());
  double nu = theta(1, 0);  // background mean

  if (isIncreasing)
    return bmd_absolute_bound(theta, BMD, BMRF * (nu - temp(0, 0)), isIncreasing);
  else
    return bmd_absolute_bound(theta, BMD, BMRF * (temp(0, 0) - nu), isIncreasing);
}

double lognormalPOLYNOMIAL_BMD_NC::bmd_point_bound(
    Eigen::MatrixXd theta, double BMD, double BMRF, bool isIncreasing
) {
  // note isIncreasing is ignored as point is specified by the user
  Eigen::MatrixXd d(1, 1);
  d(0, 0) = BMD;
  Eigen::MatrixXd temp = mean(theta, d);
  temp = exp(temp.array());

  return log(temp(0, 0)) - log(BMRF);
}

////////////////////////////////////////////////////////////////////////
// Function:  double lognormalPOLYNOMIAL_BMD_NC::bmd_hybrid_extra(Eigen::MatrixXd theta, double
// BMRF, bool isIncreasing,,double BPROB) Variables: theta - matrix of theta values for the model
//			  BMRF  - This is a value between 0 and 1 that describes the increased
// probability over BPROB
//            isIncreasing - is the function an Increasing function or decreasing function?
//            BPROB - Background probability at dose 0 considered adverse
// Purpose:   Compute the Hybrid BMD version of the hill model
//
//
//
////////////////////////////////////////////////////////////////////////
double lognormalPOLYNOMIAL_BMD_NC::bmd_hybrid_extra_bound(
    Eigen::MatrixXd theta, double BMD, double BMRF, bool isIncreasing, double TAIL_PROB
) {
  ///////////////////////////////////////////////////////////////////////////////////////
  double NOT_ADVERSE_P = 1.0 - TAIL_PROB;
  Eigen::MatrixXd d(2, 1);
  d << 0.0, BMD;
  Eigen::MatrixXd mu = mean(theta, d);  // compute the mean at background an BMD
  Eigen::MatrixXd var = variance(theta, d);

  ///////////////////////////////////////////////////////////////////////////////////////
  double ct_off =
      gsl_cdf_lognormal_Pinv(isIncreasing ? NOT_ADVERSE_P : TAIL_PROB, mu(0, 0), sqrt(var(0, 0)));
  ///////////////////////////////////////////////////////////////////////////////////////
  // cout << ct_off << endl;
  double temp;

  if (isIncreasing) {
    temp = ((1.0 - gsl_cdf_lognormal_P(ct_off, mu(1, 0), sqrt(var(1, 0)))) - TAIL_PROB) /
           NOT_ADVERSE_P;
  } else {
    temp = (gsl_cdf_lognormal_P(ct_off, mu(1, 0), sqrt(var(1, 0))) - TAIL_PROB) / NOT_ADVERSE_P;
  }

  return log(temp) - log(BMRF);
}

// Functions: lognormalPOLYNOMIAL_BMD_NC::bmd_absolute(Eigen::MatrixXd theta, double BMRF, bool
// isIncreasing)
//			  lognormalPOLYNOMIAL_BMD_NC::bmd_stdev(Eigen::MatrixXd theta, double BMRF,
// bool isIncreasing) 			  lognormalPOLYNOMIAL_BMD_NC::bmd_reldev(Eigen::MatrixXd
// theta, double BMRF, bool isIncreasing)
//            lognormalPOLYNOMIAL_BMD_NC::bmd_point(Eigen::MatrixXd theta, double BMRF, bool
//            isIncreasing)
// Purpose :  return the BMD given the parameter values theta and the BMRF. Note they are  call the
// code
//            in  lognormalPOLYNOMIAL_BMD_NC::bmd_absolute(Eigen::MatrixXd theta, double BMRF, bool
//            isIncreasing)
//
double lognormalPOLYNOMIAL_BMD_NC::bmd_absolute(
    Eigen::MatrixXd theta, double BMRF, bool isIncreasing
) {
  double min, mid, max, mu_zero;
  min = 0.0;
  max = X.maxCoeff();
  ;
  mid = (min + max) * 0.5;
  Eigen::MatrixXd d(3, 1);
  d << min, mid, max;
  Eigen::MatrixXd t_mean = mean(theta, d);
  t_mean = exp(t_mean.array());
  mu_zero = t_mean(0, 0);

  int niter = 0;

  while (fabs(t_mean(2, 0) - mu_zero) - BMRF < 0 && niter < 10) {
    niter++;
    max *= 2.0;

    d << min, mid, max;
    mid = (min + max) * 0.5;
    t_mean = mean(theta, d);
    t_mean = exp(t_mean.array());
  }

  if (niter == 10) {
    return INFINITY;  // we have gone above a reasonable value for the
                      // search for the BMD return INFINITY.
  }
  double test = fabs(t_mean(1, 0) - mu_zero) - BMRF;
  while (fabs(test) > 1e-7) {  // zero in on the BMD
    if (test > 0) {
      max = mid;
    } else {
      min = mid;
    }
    mid = (min + max) * 0.5;
    d << min, mid, max;
    t_mean = mean(theta, d);
    t_mean = exp(t_mean.array());
    test = fabs(t_mean(1, 0) - mu_zero) - BMRF;
  }
  return mid;
}

double lognormalPOLYNOMIAL_BMD_NC::bmd_stdev(
    Eigen::MatrixXd theta, double BMRF, bool isIncreasing
) {
  if (!isIncreasing) BMRF *= -1.0;
  Eigen::MatrixXd d(1, 1);
  d(0, 0) = 0.0;
  Eigen::MatrixXd temp = variance(theta, d);
  Eigen::MatrixXd med = mean(theta, d);
  med = exp(med.array());
  double t = temp(0, 0);
  Eigen::MatrixXd md = abs(exp(log(med.array()) + BMRF * pow(t, 0.5)) - med.array());
  return bmd_absolute(theta, md(0, 0), isIncreasing);
}

double lognormalPOLYNOMIAL_BMD_NC::bmd_reldev(
    Eigen::MatrixXd theta, double BMRF, bool isIncreasing
) {
  Eigen::MatrixXd d(1, 1);
  d(0, 0) = 0.0;
  Eigen::MatrixXd temp = mean(theta, d);
  temp = exp(temp.array());
  double t;

  if (isIncreasing)
    t = BMRF * temp(0, 0);
  else
    t = temp(0, 0) - BMRF * temp(0, 0);

  return bmd_absolute(theta, t, isIncreasing);
}

double lognormalPOLYNOMIAL_BMD_NC::bmd_point(
    Eigen::MatrixXd theta, double BMRF, bool isIncreasing
) {
  double min, mid, max, mu_zero;
  min = 0.0;
  max = X.maxCoeff();
  ;
  mid = (min + max) * 0.5;
  Eigen::MatrixXd d(3, 1);
  d << min, mid, max;
  Eigen::MatrixXd t_mean = mean(theta, d);
  t_mean = exp(t_mean.array());
  mu_zero = t_mean(0, 0);

  int niter = 0;

  while (isIncreasing ? t_mean(2, 0) - BMRF < 0 : t_mean(2, 0) - BMRF > 0 && niter < 10) {
    niter++;
    max *= 2.0;

    d << min, mid, max;
    mid = (min + max) * 0.5;
    t_mean = mean(theta, d);
    t_mean = exp(t_mean.array());
  }

  if (niter == 10 || isIncreasing ? t_mean(0, 0) > BMRF
                                  : t_mean(0, 0) < BMRF) {  // depends on if we are
                                                            // increasing or decreasing
    return INFINITY;  // we have gone above a reasonable value for the
                      // search for the BMD return INFINITY.
  }
  double test = t_mean(1, 0) - BMRF;
  while (fabs(test) > 1e-7) {  // zero in on the BMD
    if (isIncreasing ? test > 0 : test < 0) {
      max = mid;
    } else {
      min = mid;
    }
    mid = (min + max) * 0.5;
    d << min, mid, max;
    t_mean = mean(theta, d);
    t_mean = exp(t_mean.array());
    test = t_mean(1, 0) - BMRF;
  }
  return mid;
}

double lognormalPOLYNOMIAL_BMD_NC::bmd_extra(
    Eigen::MatrixXd theta, double BMRF, bool isIncreasing
) {
  Eigen::MatrixXd d(1, 1);
  d(0, 0) = 0.0;
  Eigen::MatrixXd temp = mean(theta, d);

  return 0.0;
}

////////////////////////////////////////////////////////////////////////
// Function:  double lognormalPOLYNOMIAL_BMD_NC::bmd_hybrid_extra(Eigen::MatrixXd theta, double
// BMRF, bool isIncreasing,,double BPROB) Variables: theta - matrix of theta values for the model
//			  BMRF  - This is a value between 0 and 1 that describes the increased
// probability over BPROB
//            isIncreasing - is the function an Increasing function or decreasing function?
//            BPROB - Background probability at dose 0 considered adverse
// Purpose:   Compute the Hybrid BMD version of the hill model
//
//
//
////////////////////////////////////////////////////////////////////////
double lognormalPOLYNOMIAL_BMD_NC::bmd_hybrid_extra(
    Eigen::MatrixXd theta, double BMRF, bool isIncreasing, double TAIL_PROB
) {
  double NOT_ADVERSE_P = 1.0 - TAIL_PROB;

  ////////////////////////////////////////////////////////////////////
  // Get the mean and variance at dose zero as well as a very high dose
  double min_d = 0.0;
  double max_d = X.maxCoeff();
  double mid = 0.5 * (min_d + max_d);
  Eigen::MatrixXd d(3, 1);
  d << min_d, mid, max_d;
  Eigen::MatrixXd temp_mean = mean(theta, d);
  Eigen::MatrixXd temp_var = variance(theta, d);
  //////////////////////////////////////////////////////////////////////
  double mu_zero = temp_mean(0, 0);
  double std_zero = sqrt(temp_var(0, 0));
  double ct_off = gsl_cdf_lognormal_Pinv(
      isIncreasing ? NOT_ADVERSE_P : TAIL_PROB, mu_zero, std_zero
  );  // CUTOFF AT DOSE = 0
  double P = TAIL_PROB + BMRF * NOT_ADVERSE_P;
  // double bmr_mult =  0.0; //gsl_cdf_lognormal_Pinv(NOT_ADVERSE_P - BMRF * (NOT_ADVERSE_P));

  double test_prob =
      isIncreasing
          ? 1.0 - gsl_cdf_lognormal_P(ct_off, temp_mean(2, 0), sqrt(temp_var(2, 0)))
          : gsl_cdf_lognormal_P(ct_off, temp_mean(2, 0), sqrt(temp_var(2, 0)));  // standardize
  double test = 0;

  int k = 0;
  while (test_prob < P && k < 10) {  // Go up to 2^10 times the maximum tested dose
                                     // if we cant find it after that we return infinity
    max_d *= 2;
    d << min_d, mid, max_d;
    temp_mean = mean(theta, d);
    temp_var = variance(theta, d);

    test_prob = isIncreasing
                    ? 1.0 - gsl_cdf_lognormal_P(ct_off, temp_mean(2, 0), sqrt(temp_var(2, 0)))
                    : gsl_cdf_lognormal_P(ct_off, temp_mean(2, 0), sqrt(temp_var(2, 0)));

    k++;
  }

  if (k == 10 || test_prob < P)  // have not been able to bound the BMD
  {
    return std::numeric_limits<double>::infinity();
  }

  test_prob = isIncreasing
                  ? 1.0 - gsl_cdf_lognormal_P(ct_off, temp_mean(1, 0), sqrt(temp_var(1, 0)))
                  : gsl_cdf_lognormal_P(ct_off, temp_mean(1, 0), sqrt(temp_var(1, 0)));

  double temp_test = test_prob - P;
  /////////////////////////////////////////////////////////////////////////////
  while (fabs(temp_test) > 1e-5) {
    // we have bounded the BMD now we use a root finding algorithm to
    // figure out what it is default difference is a probability of of 1e-5
    if (temp_test > 0) {
      max_d = mid;
    } else {
      min_d = mid;
    }

    mid = 0.5 * (max_d + min_d);
    d << min_d, mid, max_d;

    temp_mean = mean(theta, d);
    temp_var = variance(theta, d);

    test_prob = isIncreasing
                    ? 1.0 - gsl_cdf_lognormal_P(ct_off, temp_mean(1, 0), sqrt(temp_var(1, 0)))
                    : gsl_cdf_lognormal_P(ct_off, temp_mean(1, 0), sqrt(temp_var(1, 0)));

    temp_test = test_prob - P;
  }

  return mid;
}
