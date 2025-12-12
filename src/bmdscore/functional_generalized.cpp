// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "functional_generalized.h"

#include <random>
// we only include RcppEigen.h which pulls Rcpp.h in for us
// #include <RcppEigen.h>
// #include <math.h>
// #include <cmath>
// #include <boost/math/distributions/gamma.hpp>
// #include <gsl/gsl_cdf.h>
// #include <gsl/gsl_randist.h>
// #include <gsl/gsl_rng.h>
//
//
// #include <Eigen/Dense>

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]

// simple example of creating two matrices and
// returning the result of an operation on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
// using namespace Rcpp;

//////////////////////////////////////////////////////////////////////////////
/// @brief Computes the  transformation
/// @param params - A (nparams x 1) vector of Q. Linear parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the Q. Linear transformed doses
Eigen::VectorXd binomial_cqlinear(const Eigen::VectorXd& params, const Eigen::MatrixXd& X) {
  double a1 = params[0];
  double b1 = params[1];
  double c1 = params[2];
  double sum = a1 + b1 + c1;
  double p_zero = a1 / sum;
  double p_one = (a1 + b1) / sum;

  double g = p_zero;  // R::qnorm(p_zero,0,1,true,false);
  double b = -log((1 - p_one) / (1 - p_zero));

  Eigen::VectorXd out = X.col(0);

  out = out.unaryExpr([&g, &b](double xx) { return g + (1 - g) * (1 - exp(-b * xx)); });
  return out;
}

/// @brief Computes the  transformation
/// @param params - A (nparams x 1) vector of multistage2 parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the multistage2 transformed doses
Eigen::VectorXd binomial_cmstage2(const Eigen::VectorXd& params, const Eigen::MatrixXd& X) {
  double a11 = params[0];
  double b11 = params[1];
  double c11 = params[2];
  double sum = a11 + b11 + c11;
  double p_zero = a11 / sum;
  double p_one = (a11 + b11) / sum;

  double y1 = params[3];

  double g = p_zero;  // R::qnorm(p_zero,0,1,true,false);
  double b1 = -(y1)*log((1 - p_one) / (1 - p_zero));
  double b2 = -(1.0 - y1) * log((1 - p_one) / (1 - p_zero));

  Eigen::VectorXd out = X.col(0);

  out = out.unaryExpr([&g, &b1, &b2](double xx) {
    return g + (1 - g) * (1 - exp(-b1 * xx - b2 * xx * xx));
  });
  return out;
}

//////////////////////////////////////////////////////////////////////////////
/// @brief Computes the  transformation
/// @param params - A (nparams x 1) vector of logprobit parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the logprobit transformed doses
Eigen::VectorXd binomial_clogprobit(const Eigen::VectorXd& params, const Eigen::MatrixXd& X) {
  double a1 = params[0];
  double b1 = params[1];
  double c1 = params[2];
  double b = params[3];
  double sum = a1 + b1 + c1;
  double p0 = a1 / sum;
  double p1 = (a1 + b1) / sum;

  //  double a = R::qnorm((p1-p0)/(1-p0),0,1,true,false);
  double a = gsl_cdf_gaussian_Pinv((p1 - p0) / (1 - p0), 1);

  Eigen::VectorXd out = X.col(0);

  out = out.unaryExpr([&p0, &b, &a](double xx) {
    if (xx == 0) {
      return p0;
    }
    // return p0+(1-p0)*R::pnorm(a+b*log(xx),0,1,true,false);
    return p0 + (1 - p0) * gsl_cdf_gaussian_P(a + b * log(xx), 1);
  });
  return out;
}

/// @brief
/// @param params - A (nparams x 1) vector of loglogistic parameters
/// @param X
/// @return
Eigen::VectorXd binomial_cloglogistic(const Eigen::VectorXd& params, const Eigen::MatrixXd& X) {
  double a1 = params[0];
  double b1 = params[1];
  double c1 = params[2];
  double b = params[3];
  double sum = a1 + b1 + c1;
  double p0 = a1 / sum;
  double p1 = (a1 + b1) / sum;
  double a = -log((1 - p1) / (p1 - p0));

  Eigen::VectorXd out = X.col(0);
  if (b <= 0) {
    b = 1e-8;
  }

  out = out.unaryExpr([&p0, &b, &a](double xx) {
    if (xx == 0) {
      return p0;
    }
    return p0 + (1 - p0) / (1 + exp(-a - b * log(xx)));
  });
  return out;
}

/// @brief Computes the  transformation
/// @param params - A (nparams x 1) vector of Probit parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the probit transformed doses
Eigen::VectorXd binomial_cprobit(const Eigen::VectorXd& params, const Eigen::MatrixXd& X) {
  double a1 = params[0];
  double b1 = params[1];
  double c1 = params[2];
  double sum = a1 + b1 + c1;
  double p_zero = a1 / sum;
  double p_one = (a1 + b1) / sum;

  //  double a = R::qnorm(p_zero,0,1,true,false);
  //  double b = R::qnorm(p_one ,0,1,true,false) - a;
  double a = gsl_cdf_gaussian_Pinv(p_zero, 1);
  double b = gsl_cdf_gaussian_Pinv(p_one, 1) - a;

  Eigen::VectorXd out = X.col(0);

  // out = out.unaryExpr([&a,&b](double xx){return R::pnorm(a+b*xx,0,1,true,false);});
  out = out.unaryExpr([&a, &b](double xx) { return gsl_cdf_gaussian_P(a + b * xx, 1); });
  return out;
}

/// @brief Computes the  transformation
/// @param params - A (nparams x 1) vector of D. Hill parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the D. Hill transformed doses
Eigen::VectorXd binomial_chill(const Eigen::VectorXd& params, const Eigen::MatrixXd& X) {
  double g = params[0];
  double v = params[1];
  double a = params[2];
  double b = params[3];

  Eigen::VectorXd out = X.col(0);

  out = out.unaryExpr([&g, &b, &a, &v](double xx) {
    if (xx == 0) {
      return g;
    }
    return g + v * (1 - g) / (1 + exp(-a - b * log(xx)));
  });

  // Rcout << out << '\n';
  return out;
}

/// @brief The Weibull Model for dichotomous data
/// @param params - A (nparams x 1) vector of Weibull parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of probability
Eigen::VectorXd binomial_cweibull(const Eigen::VectorXd& params, const Eigen::MatrixXd& X) {
  double a1 = params[0];
  double b1 = params[1];
  double c1 = params[2];
  double sum = a1 + b1 + c1;
  double p_zero = a1 / sum;
  double p_one = (a1 + b1) / sum;
  double g = params[3];

  double a = p_zero;
  double b = -log((1 - p_one) / (1 - p_zero));

  Eigen::VectorXd out = X.col(0);

  out = out.unaryExpr([&a, &b, &g](double xx) { return a + (1 - a) * (1 - exp(-b * pow(xx, g))); });

  return out;
}

/// @brief The Logistic Model for dichotomous data
/// @param params - A (nparams x 1) vector of Logistic parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of probability
Eigen::VectorXd binomial_clogistic(const Eigen::VectorXd& params, const Eigen::MatrixXd& X) {
  double a1 = params[0];
  double b1 = params[1];
  double c1 = params[2];
  double sum = a1 + b1 + c1;
  double p_zero = a1 / sum;
  double p_one = (a1 + b1) / sum;

  double a = log(p_zero / (1 - p_zero));
  double b = log(p_one / (1 - p_one)) - a;

  Eigen::VectorXd out = X.col(0);

  out = out.unaryExpr([&a, &b](double xx) { return 1 / (1 + exp(-a - b * xx)); });
  return out;
}

/// @brief The Gamma Model for dichotomous data
/// @param params - A (nparams x 1) vector of Gamma parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of probability
Eigen::VectorXd binomial_cgamma(const Eigen::VectorXd& params, const Eigen::MatrixXd& X) {
  double a1 = params[0];
  double b1 = params[1];
  double c1 = params[2];
  double sum = a1 + b1 + c1;
  double p0 = a1 / sum;
  double p1 = (a1 + b1) / sum;
  double alpha = params[3];

  if (alpha <= 0.2) {
    alpha = 0.2;
  }

  // double b = R::qgamma((p1-p0)/(1-p0),alpha,1,true, false);
  double b = gsl_cdf_gamma_Pinv((p1 - p0) / (1 - p0), alpha, 1);

  Eigen::VectorXd out = X.col(0);
  if (b <= 0) {
    b = 1e-8;
  }
  // out = out.unaryExpr([&p0,&b,&alpha](double xx){return p0+(1-p0)*(
  // R::pgamma(b*xx,alpha,1.0,true,false));});
  out = out.unaryExpr([&p0, &b, &alpha](double xx) {
    return p0 + (1 - p0) * (gsl_cdf_gamma_P(b * xx, alpha, 1.0));
  });

  return out;
}

/// @brief
/// @param x
/// @param alpha
/// @param beta
/// @return
double log_beta_pdf(double x, double alpha, double beta) {
  double returnV = 0.0;
  if (x <= 0 || x >= 1.0) {
    return -std::numeric_limits<double>::infinity();
  }
  returnV += (alpha - 1.0) * log(x);
  returnV += (beta - 1.0) * log(1.0 - x);

  return returnV;
}

double log_gamma_pdf(double x, double alpha, double beta) {
  double returnV = 0.0;

  if (x <= 0) {
    return -std::numeric_limits<double>::infinity();
  }
  returnV += (alpha - 1.0) * log(x);
  returnV += -(beta)*x;

  return returnV;
}

/// @brief Computes the  transformation
/// @param params - A (nparams x 1) vector of multistage3 parameters (DO NOT USE)
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the multistage3 transformed doses
Eigen::VectorXd binomial_mstage3(const Eigen::VectorXd& params, const Eigen::MatrixXd& X) {
  double p_zero = params[0];
  double p_one = params[1];
  double a1 = params[2];
  double a2 = params[3];
  double a3 = params[4];

  double sum_a = a1 + a2 + a3;
  double d1 = a1 / sum_a;
  double d2 = a2 / sum_a;
  double d3 = 1 - (d1 + d2);

  double g = p_zero;  // R::qnorm(p_zero,0,1,true,false);
  double b1 = -d1 * log((1 - p_one) / (1 - p_zero));
  double b2 = -d2 * log((1 - p_one) / (1 - p_zero));
  double b3 = -d3 * log((1 - p_one) / (1 - p_zero));

  Eigen::VectorXd out = X.col(0);

  out = out.unaryExpr([&g, &b1, &b2, &b3](double xx) {
    return g + (1 - g) * (1 - exp(-b1 * xx - b2 * xx * xx - b3 * xx * xx * xx));
  });
  return out;
}

/// @brief Computes the Hill transformation (DO NOT USE)
/// @param params - A (nparams x 1) vector of Hill parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the Hill transformed doses
Eigen::VectorXd hill_cpp2(const Eigen::VectorXd& params, const Eigen::MatrixXd& X) {
  double a = params[0];
  double b = params[1];
  double c = params[2];
  double d = params[3];
  Eigen::VectorXd out = X.col(0);
  out = out.unaryExpr([&a, &b, &c, &d](double xx) {
    return a * (1.0 + (c - 1.0) * (1.0 - pow(b, d) / (pow(b, d) + pow(xx, d))));
  });
  return out;
}

/// @brief Computes the Power transformation (DO NOT USE)
/// @param params - A (nparams x 1) vector of Power parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the Power transformed doses
Eigen::VectorXd power_cpp2(const Eigen::VectorXd& params, const Eigen::MatrixXd& X) {
  double a = params[0];
  double b = params[1];
  double c = params[2];
  Eigen::VectorXd out = a + b * pow(X.col(0).array(), c);
  return out;
}

/// @brief Computes the Exponential transformation (DO NOT USE)
/// @param params - A (nparams x 1) vector of Exponential parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the Exponential transformed doses
Eigen::VectorXd exponential_cpp2(const Eigen::VectorXd& params, const Eigen::MatrixXd& X) {
  double a = params[0];
  double b = params[1];
  double c = params[2];
  double d = params[3];
  Eigen::VectorXd out = X.col(0);
  out = out.unaryExpr([&a, &b, &c, &d](double xx) {
    return a * (1.0 + (c - 1.0) * (1.0 - exp(-b * pow(xx, d))));
  });
  return out;
}

/// @brief Computes the Linear transformation (DO NOT USE)
/// @param params - A (nparams x 1) vector of Linear parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the Linear transformed doses
Eigen::VectorXd linear_cpp2(const Eigen::VectorXd& params, const Eigen::MatrixXd& X) {
  double a = params[0];
  double b = params[1];
  Eigen::VectorXd out = a + b * X.col(0).array();
  return out;
}

/// @brief Computes the Linear transformation
/// @param params - A (nparams x 1) vector of Linear parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the Linear transformed doses
Eigen::VectorXd continuous_linear_transform(
    const Eigen::VectorXd& params, const Eigen::MatrixXd& X
) {
  double m0 = params[0];  // predicted mean at dose = 0
  double m1 = params[1];  // predicted mean at dose = 1

  double g = m0;
  double b = m1 - m0;
  Eigen::VectorXd out = X.col(0);
  out = out.unaryExpr([&g, &b](double xx) { return g + b * xx; });
  return out;
}

/// @brief Computes the Poly2 transformation
/// @param params - A (nparams x 1) vector of Poly2 parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the Poly2 transformed doses
Eigen::VectorXd continuous_poly2_transform(
    const Eigen::VectorXd& params, const Eigen::MatrixXd& X
) {
  double m0 = params[0];  // mean at dose 0
  double m1 = params[1];  // mean at dose 1
  double y1 = params[2];  // unconstrained split param (matches prior_v now)

  // Direction and magnitude of total change
  double delta = m1 - m0;
  double dir = (delta >= 0.0) ? 1.0 : -1.0;
  double mag = std::fabs(delta);

  double b1 = dir * y1 * mag;
  double b2 = dir * (1.0 - y1) * mag;

  // Defensive sign enforcement (should already hold, but keeps things robust)
  if (dir > 0.0) {  // increasing
    if (b1 <= 0.0) b1 = 1e-08;
    if (b2 <= 0.0) b2 = 1e-08;
  } else {  // decreasing
    if (b1 >= 0.0) b1 = -1e-08;
    if (b2 >= 0.0) b2 = -1e-08;
  }

  // Evaluate f(x) = m0 + b1 x + b2 x^2
  const double g = m0;
  Eigen::VectorXd out = X.col(0);
  out = out.unaryExpr([&](double xx) { return g + b1 * xx + b2 * xx * xx; });
  return out;
}

/// @brief Computes the Power transformation
/// @param params - A (nparams x 1) vector of Power parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the Power transformed doses
Eigen::VectorXd continuous_power_transform(
    const Eigen::VectorXd& params, const Eigen::MatrixXd& X
) {
  double m0 = params[0];  // response at dose = 0
  double m1 = params[1];  // response at dose = 1
  double n = params[2];   // power parameter

  double g = m0;
  double v = m1 - m0;

  Eigen::VectorXd out = X.col(0);
  out = out.unaryExpr([&g, &v, &n](double xx) { return g + v * std::pow(xx, n); });

  return out;
}

/// @brief Computes the Exp3 transformation
/// @param params - A (nparams x 1) vector of Exp3 parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the Exp3 transformed doses
Eigen::VectorXd continuous_exp3_transform(const Eigen::VectorXd& params, const Eigen::MatrixXd& X) {
  double m0 = params[0];  // response at dose = 0
  double m1 = params[1];  // response at dose = 1
  double c = params[2];   // power shape

  double log_ratio = std::log(m1 / m0);
  double sign = (log_ratio >= 0) ? 1.0 : -1.0;
  double b = std::pow(std::abs(log_ratio), 1.0 / c);

  Eigen::VectorXd out = X.col(0);
  out = out.unaryExpr([&m0, &b, &c, &sign](double xx) {
    return m0 * std::exp(sign * std::pow(b * xx, c));
  });

  return out;
}

/// @brief Computes the Exp3-log transformation
/// @param params - A (nparams x 1) vector of Exp3-log parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the Exp3-log transformed doses
Eigen::VectorXd continuous_exp3_transform_log(
    const Eigen::VectorXd& params, const Eigen::MatrixXd& X
) {
  double m0 = std::exp(params[0]);  // response at dose = 0
  double m1 = std::exp(params[1]);  // response at dose = 1
  double c = params[2];             // power shape

  double log_ratio = std::log(m1 / m0);
  double sign = (log_ratio >= 0) ? 1.0 : -1.0;
  double b = std::pow(std::abs(log_ratio), 1.0 / c);

  Eigen::VectorXd out = X.col(0);
  out = out.unaryExpr([&m0, &b, &c, &sign](double xx) {
    return m0 * std::exp(sign * std::pow(b * xx, c));
  });

  return out.array().log().matrix();
}

/// @brief Computes the Exp5 transformation
/// @param params - A (nparams x 1) vector of Exp5 parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the Exp5 transformed doses
Eigen::VectorXd continuous_exp5_transform(const Eigen::VectorXd& params, const Eigen::MatrixXd& X) {
  double m0 = params[0];  // response at dose = 0
  double m1 = params[1];  // response at dose = 1
  double b = params[2];   // slope-like parameter
  double n = params[3];   // power parameter

  double exp_bn = std::exp(std::pow(b, n));
  double numerator = m0 - m1 * exp_bn;
  double denominator = m0 - m0 * exp_bn;

  double c = numerator / denominator;
  if (c <= 0.0 || !std::isfinite(c)) {
    c = 1e-8;
  }

  Eigen::VectorXd out = X.col(0);
  out = out.unaryExpr([&m0, &b, &n, &c](double xx) {
    return m0 * (c - (c - 1.0) * std::exp(-std::pow(b * xx, n)));
  });

  return out;
}

/// @brief Computes the Exp5-log transformation
/// @param params - A (nparams x 1) vector of Exp5-log parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the Exp5-log transformed doses
Eigen::VectorXd continuous_exp5_transform_log(
    const Eigen::VectorXd& params, const Eigen::MatrixXd& X
) {
  double m0 = std::exp(params[0]);  // response at dose = 0
  double m1 = std::exp(params[1]);  // response at dose = 1
  double b = params[2];             // slope-like parameter
  double n = params[3];             // power parameter

  double exp_bn = std::exp(std::pow(b, n));
  double numerator = m0 - m1 * exp_bn;
  double denominator = m0 - m0 * exp_bn;

  double c = numerator / denominator;
  if (c <= 0.0 || !std::isfinite(c)) {
    c = 1e-8;
  }

  Eigen::VectorXd out = X.col(0);
  out = out.unaryExpr([&m0, &b, &n, &c](double xx) {
    return m0 * (c - (c - 1.0) * std::exp(-std::pow(b * xx, n)));
  });

  return out.array().log().matrix();
}

/// @brief Computes the BMDS Hill transformation
/// @param params - A (nparams x 1) vector of BMDS Hill parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the BMDS Hill transformed doses
Eigen::VectorXd continuous_hill_transform(const Eigen::VectorXd& params, const Eigen::MatrixXd& X) {
  double m0 = params[0];  // response at dose = 0
  double m1 = params[1];  // response at dose = 1
  double k = params[2];   // half-saturation constant
  double n = params[3];   // Hill coefficient (power)

  double kn = std::pow(k, n);
  double v = (m1 - m0) * (kn + 1);

  Eigen::VectorXd out = X.col(0);
  out = out.unaryExpr([&m0, &v, &k, &n](double d) {
    double dn = std::pow(d, n);
    return m0 + v * dn / (std::pow(k, n) + dn);
  });

  return out;
}

/// @brief Computes the EFSA Hill transformation
/// @param params - A (nparams x 1) vector of EFSA Hill parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the EFSA Hill transformed doses
Eigen::VectorXd continuous_hill4_efsa_transform(
    const Eigen::VectorXd& params, const Eigen::MatrixXd& X
) {
  double m0 = params[0];  // response at dose = 0
  double m1 = params[1];  // response at dose = 1
  double b = params[2];   // half-saturation constant
  double d = params[3];   // Hill coefficient (power)

  double bd = std::pow(b, d);
  double c = ((m1 - m0) / m0) * (bd + 1.0) + 1.0;

  Eigen::VectorXd out = X.col(0);
  out = out.unaryExpr([&m0, &c, &bd, &d](double x) {
    double xd = std::pow(x, d);
    return m0 * (1.0 + (c - 1.0) * xd / (bd + xd));
  });

  return out;
}

/// @brief Computes the EFSA Hill-log transformation
/// @param params - A (nparams x 1) vector of EFSA Hill-log parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the EFSA Hill-log transformed doses
Eigen::VectorXd continuous_hill4_efsa_transform_log(
    const Eigen::VectorXd& params, const Eigen::MatrixXd& X
) {
  double m0 = std::exp(params[0]);
  double m1 = std::exp(params[1]);
  double b = params[2];
  double d = params[3];

  double bd = std::pow(b, d);
  double c = ((m1 - m0) / m0) * (bd + 1.0) + 1.0;

  Eigen::VectorXd out = X.col(0);
  out = out.unaryExpr([&m0, &c, &bd, &d](double x) {
    double xd = std::pow(x, d);
    return m0 * (1.0 + (c - 1.0) * xd / (bd + xd));
  });

  return out.array().log().matrix();
}

/// @brief Computes the EFSA Inverse Exponential transformation
/// @param params - A (nparams x 1) vector of EFSA Inverse Exponential parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the EFSA Inverse Exponential transformed doses
Eigen::VectorXd continuous_invexp_efsa_transform(
    const Eigen::VectorXd& params, const Eigen::MatrixXd& X
) {
  double m0 = params[0];  // response at dose = 0
  double m1 = params[1];  // response at dose = 1
  double b = params[2];   // shape parameter
  double d = params[3];   // power parameter

  double exp_neg_b = std::exp(-b);
  double c = (m1 - m0 + m0 * exp_neg_b) / (m0 * exp_neg_b);

  Eigen::VectorXd out = X.col(0);
  out = out.unaryExpr([&m0, &c, &b, &d](double x) {
    double xd = std::pow(x, -d);
    return m0 * (1.0 + (c - 1.0) * (std::exp(-b * xd)));
  });

  return out;
}

/// @brief Computes the EFSA Inverse Exponential-log transformation
/// @param params - A (nparams x 1) vector of EFSA Inverse Exponential-log parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the EFSA Inverse Exponential-log transformed doses
Eigen::VectorXd continuous_invexp_efsa_transform_log(
    const Eigen::VectorXd& params, const Eigen::MatrixXd& X
) {
  double m0 = std::exp(params[0]);
  double m1 = std::exp(params[1]);
  double b = params[2];
  double d = params[3];

  double exp_neg_b = std::exp(-b);
  double c = (m1 - m0 + m0 * exp_neg_b) / (m0 * exp_neg_b);

  Eigen::VectorXd out = X.col(0);
  out = out.unaryExpr([&m0, &c, &b, &d](double x) {
    double xd = std::pow(x, -d);
    return m0 * (1.0 + (c - 1.0) * (std::exp(-b * xd)));
  });

  return out.array().log().matrix();
}

/// @brief Computes the EFSA Lognormal transformation
/// @param params - A (nparams x 1) vector of EFSA Lognormal parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the EFSA Lognormal transformed doses
Eigen::VectorXd continuous_lognormal_efsa_transform(
    const Eigen::VectorXd& params, const Eigen::MatrixXd& X
) {
  double m0 = params[0];
  double m1 = params[1];
  double b = params[2];
  double d = params[3];

  double logb = std::log(b);
  // double phi_logb = R::pnorm(logb, 0.0, 1.0, true, false);
  double phi_logb = gsl_cdf_gaussian_P(logb, 1.0);

  double c = (m1 - m0 + m0 * phi_logb) / (m0 * phi_logb);

  Eigen::VectorXd out = X.col(0);
  out = out.unaryExpr([&m0, &c, &b, &d](double x) {
    if (x == 0) {
      return m0;
    }
    double z = std::log(b) + d * std::log(x);
    // double phi = R::pnorm(z, 0.0, 1.0, true, false);
    double phi = gsl_cdf_gaussian_P(z, 1.0);
    return m0 * (1.0 + (c - 1.0) * phi);
  });

  return out;
}

/// @brief Computes the EFSA Lognormal-log transformation
/// @param params - A (nparams x 1) vector of EFSA Lognormal-log parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the EFSA Lognormal-log transformed doses
Eigen::VectorXd continuous_lognormal_efsa_transform_log(
    const Eigen::VectorXd& params, const Eigen::MatrixXd& X
) {
  double m0 = std::exp(params[0]);
  double m1 = std::exp(params[1]);
  double b = params[2];
  double d = params[3];

  // double phi_logb = R::pnorm(std::log(b), 0.0, 1.0, true, false);
  double phi_logb = gsl_cdf_gaussian_P(std::log(b), 1.0);
  double c = (m1 - m0 + m0 * phi_logb) / (m0 * phi_logb);

  Eigen::VectorXd out = X.col(0);
  out = out.unaryExpr([&m0, &c, &b, &d](double x) {
    if (x == 0) return m0;
    double z = std::log(b) + d * std::log(x);
    // double phi = R::pnorm(z, 0.0, 1.0, true, false);
    double phi = gsl_cdf_gaussian_P(z, 1.0);
    return m0 * (1.0 + (c - 1.0) * phi);
  });

  return out.array().log().matrix();
}

/// @brief Computes the EFSA Gamma transformation
/// @param params - A (nparams x 1) vector of EFSA Gamma parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the EFSA Gamma transformed doses
Eigen::VectorXd continuous_gamma_efsa_transform(
    const Eigen::VectorXd& params, const Eigen::MatrixXd& X
) {
  double m0 = params[0];
  double m1 = params[1];
  double b = params[2];
  double d = params[3];

  // double pg = R::pgamma(b, d, 1.0, true, false);  // pgamma(b, d, scale=1)
  double pg = gsl_cdf_gamma_P(b, d, 1.0);  // pgamma(b, d, scale=1)
  double c = (m1 - m0) / (m0 * pg) + 1.0;

  Eigen::VectorXd out = X.col(0);
  out = out.unaryExpr([&m0, &c, &b, &d](double x) {
    double bx = b * x;
    // double pg = R::pgamma(bx, d, 1.0, true, false);
    double pg = gsl_cdf_gamma_P(bx, d, 1.0);
    return m0 * (1.0 + (c - 1.0) * pg);
  });

  return out;
}

/// @brief Computes the EFSA Gamma-log transformation
/// @param params - A (nparams x 1) vector of EFSA Gamma-log parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the EFSA Gamma-log transformed doses
Eigen::VectorXd continuous_gamma_efsa_transform_log(
    const Eigen::VectorXd& params, const Eigen::MatrixXd& X
) {
  double m0 = std::exp(params[0]);
  double m1 = std::exp(params[1]);
  double b = params[2];
  double d = params[3];

  // double pg = R::pgamma(b, d, 1.0, true, false);
  double pg = gsl_cdf_gamma_P(b, d, 1.0);
  double c = (m1 - m0) / (m0 * pg) + 1.0;

  Eigen::VectorXd out = X.col(0);
  out = out.unaryExpr([&m0, &c, &b, &d](double x) {
    double bx = b * x;
    // double pg = R::pgamma(bx, d, 1.0, true, false);
    double pg = gsl_cdf_gamma_P(bx, d, 1.0);
    return m0 * (1.0 + (c - 1.0) * pg);
  });

  return out.array().log().matrix();
}

/// @brief Computes the EFSA LMS2 transformation
/// @param params - A (nparams x 1) vector of EFSA LMS2 parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the EFSA LMS2 transformed doses
Eigen::VectorXd continuous_lms_efsa_transform(
    const Eigen::VectorXd& params, const Eigen::MatrixXd& X
) {
  double m0 = params[0];  // a
  double m1 = params[1];
  double b = params[2];
  double d = params[3];

  double exp_neg_bd = std::exp(-b - d);
  double c = (m1 - m0 * exp_neg_bd) / (m0 * (1.0 - exp_neg_bd));

  Eigen::VectorXd out = X.col(0);
  out = out.unaryExpr([&m0, &c, &b, &d](double x) {
    double exp_term = std::exp(-b * x - d * x * x);
    return m0 * (1.0 + (c - 1.0) * (1.0 - exp_term));
  });

  return out;
}

/// @brief Computes the EFSA LMS2-log transformation
/// @param params - A (nparams x 1) vector of EFSA LMS2-log parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the EFSA LMS2-log transformed doses
Eigen::VectorXd continuous_lms_efsa_transform_log(
    const Eigen::VectorXd& params, const Eigen::MatrixXd& X
) {
  double m0 = std::exp(params[0]);
  double m1 = std::exp(params[1]);
  double b = params[2];
  double d = params[3];

  double exp_neg_bd = std::exp(-b - d);
  double c = (m1 - m0 * exp_neg_bd) / (m0 * (1.0 - exp_neg_bd));

  Eigen::VectorXd out = X.col(0);
  out = out.unaryExpr([&m0, &c, &b, &d](double x) {
    double exp_term = std::exp(-b * x - d * x * x);
    return m0 * (1.0 + (c - 1.0) * (1.0 - exp_term));
  });

  return out.array().log().matrix();
}

/// @brief Computes the Logistic transformation (DO NOT USE)
/// @param params - A (nparams x 1) vector of Logistic parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the Logistic transformed doses
Eigen::VectorXd logistic_cpp2(const Eigen::VectorXd& params, const Eigen::MatrixXd& X) {
  double a = params[0];
  double b = params[1];
  double c = params[2];
  double d = params[3];
  Eigen::VectorXd out = X.col(0);
  out = out.unaryExpr([&a, &b, &c, &d](double xx) {
    return c / (1.0 + exp(-1.0 * a - b * pow(xx, d)));
  });
  return out;
}

/// @brief Creates the pointer to the correct nonlinearity
/// @param op - An integer; 1 = Hill, 2 = Power, 3 = Exponential, 5 = Logistic, otherwise use Linear
/// @return The function pointer
ptr2 choose_nonlinearity2(int op) {
  if (op == 1) {
    return &hill_cpp2;
  } else if (op == 2) {
    return &power_cpp2;
  } else if (op == 3) {
    return &exponential_cpp2;
  } else if (op == 5) {
    return &logistic_cpp2;
  } else if (op == 10) {
    return &binomial_clogistic;
  } else if (op == 11) {
    return &binomial_cweibull;
  } else if (op == 12) {
    return &binomial_cgamma;
  } else if (op == 13) {
    return &binomial_chill;
  } else if (op == 14) {
    return &binomial_cloglogistic;
  } else if (op == 15) {
    return &binomial_clogprobit;
  } else if (op == 16) {
    return &binomial_cprobit;
  } else if (op == 17) {
    return &binomial_cqlinear;
  } else if (op == 18) {
    return &binomial_cmstage2;
  } else if (op == 19) {
    return &binomial_mstage3;

    // Reparameterized continuous models
  } else if (op == 101) {
    return &continuous_linear_transform;
  } else if (op == 102) {
    return &continuous_exp3_transform;
  } else if (op == 107) {
    return &continuous_exp3_transform_log;
  } else if (op == 103) {
    return &continuous_power_transform;
  } else if (op == 104) {
    return &continuous_exp5_transform;
  } else if (op == 105) {
    return &continuous_exp5_transform_log;
  } else if (op == 106) {
    return &continuous_hill_transform;
  } else if (op == 108) {
    return &continuous_hill4_efsa_transform;
  } else if (op == 109) {
    return &continuous_hill4_efsa_transform_log;
  } else if (op == 110) {
    return &continuous_invexp_efsa_transform;
  } else if (op == 111) {
    return &continuous_invexp_efsa_transform_log;
  } else if (op == 112) {
    return &continuous_lognormal_efsa_transform;
  } else if (op == 113) {
    return &continuous_lognormal_efsa_transform_log;
  } else if (op == 114) {
    return &continuous_gamma_efsa_transform;
  } else if (op == 115) {
    return &continuous_gamma_efsa_transform_log;
  } else if (op == 116) {
    return &continuous_lms_efsa_transform;
  } else if (op == 117) {
    return &continuous_lms_efsa_transform_log;
  } else if (op == 118) {
    return &continuous_poly2_transform;

  } else {
    return &linear_cpp2;
  }
}

// from Dirk, https://github.com/RcppCore/Rcpp/issues/967
// Rcpp::NumericVector Quantile3(Rcpp::NumericVector x,
//                              Rcpp::NumericVector probs) {
Eigen::VectorXd Quantile3(Eigen::VectorXd x, Eigen::VectorXd probs) {
  // implementation of type 7
  const size_t n = x.size(), np = probs.size();
  if (n == 0) return x;
  if (np == 0) return probs;
  ////  Rcpp::NumericVector index = (n-1.)*probs, y=x.sort(), x_hi(np), qs(np);
  ////  Rcpp::NumericVector lo = Rcpp::floor(index), hi = Rcpp::ceiling(index);
  //  std::vector<double> index(probs.size());
  //  std::vector<int> lo(probs.size());
  //  std::vector<int> hi(probs.size());
  //  for (int i=0; i<probs.size(); i++){
  //    index[i] = (n-1.)*probs[i];
  //    lo[i] = floor(index[i]);
  //    hi[i] = ceil(index[i]);
  //  }

  Eigen::VectorXd index = (n - 1.) * probs;
  Eigen::VectorXd lo = Eigen::floor(index.array());
  Eigen::VectorXd hi = Eigen::ceil(index.array());

  //  std::vector<double> y= x;
  //  std::sort(y.begin(), y.end());
  //  std::vector<double> x_hi(np);
  //  std::vector<double> qs(np);
  Eigen::VectorXd y = x;
  std::sort(y.begin(), y.end());
  Eigen::VectorXd x_hi(np);
  Eigen::VectorXd qs(np);

  for (size_t i = 0; i < np; ++i) {
    qs[i] = y[lo[i]];
    x_hi[i] = y[hi[i]];
    if ((index[i] > lo[i]) && (x_hi[i] != qs[i])) {
      double h;
      h = index[i] - lo[i];
      qs[i] = (1. - h) * qs[i] + h * x_hi[i];
    }
  }
  return qs;
}

Eigen::MatrixXd truncated_linear_cpp3(const Eigen::VectorXd& x, const Eigen::VectorXd& knots) {
  Eigen::MatrixXd m(x.size(), knots.size() + 2);
  m.col(0) = Eigen::VectorXd::Ones(x.size());
  m.col(1) = x;
  for (auto j = 0; j < knots.size(); j++) {
    Eigen::VectorXd temp = x.array() - knots[j];
    m.col(j + 2) = temp.unaryExpr([](double ll) { return ll > 0.0 ? ll * ll : 0.0; });
  }
  return m;
}

Eigen::MatrixXd trunc_lin_numeric3(const double& x, const Eigen::VectorXd knots) {
  Eigen::MatrixXd m(1, knots.size() + 2);
  m(0, 0) = 1.0;
  m(0, 1) = x;
  for (auto j = 0; j < knots.size(); j++) {
    double temp = x - knots[j];
    m(0, j + 2) = temp > 0.0 ? temp * temp : 0.0;
  }

  return m;
}

// List compute_transform_f_lag1_cpp3(const Eigen::MatrixXd Y,
//                                    const Rcpp::NumericVector qtiles) {
void compute_transform_f_lag1_cpp3(
    const Eigen::MatrixXd Y, const Eigen::VectorXd qtiles, Eigen::MatrixXd& beta_return,
    Eigen::MatrixXd& knot_return
) {
  //  Eigen::MatrixXd beta_return = Eigen::MatrixXd::Zero(qtiles.size() + 2, Y.cols());
  //  Eigen::MatrixXd knot_return = Eigen::MatrixXd::Zero(qtiles.size(), Y.cols());
  beta_return = Eigen::MatrixXd::Zero(qtiles.size() + 2, Y.cols());
  knot_return = Eigen::MatrixXd::Zero(qtiles.size(), Y.cols());
  // beta_return.setZero();
  // knot_return.setZero();

  beta_return(0, 0) = Y.col(0).mean();

  for (auto jj = 1; jj < Y.cols(); jj++) {
    Eigen::VectorXd temp = Y.col(jj - 1);
    // Eigen::VectorXd knots = <Eigen::VectorXd>(Quantile3(wrap(temp), qtiles));
    Eigen::VectorXd knots = Quantile3(temp, qtiles);

    Eigen::MatrixXd tx = truncated_linear_cpp3(Y.col(jj - 1), knots);

    Eigen::MatrixXd tempp = tx.transpose() * tx;
    // add fudge factor to the diagonal
    tempp.diagonal().array() += 0.001;
    beta_return.col(jj) = tempp.colPivHouseholderQr().solve(tx.transpose() * Y.col(jj));
    knot_return.col(jj) = knots;
  }

  //  return List::create(Named("betas") = beta_return,
  //                      Named("knots") = knot_return);
  return;
}

// List compute_cov_eta_cpp3(Eigen::MatrixXd Y,
//                           List beta_return) {
void compute_cov_eta_cpp3(
    Eigen::MatrixXd Y, Eigen::MatrixXd& betas, Eigen::MatrixXd& knots, Eigen::VectorXd& colm,
    Eigen::MatrixXd& tempp
) {
  Eigen::MatrixXd eta = Y;

  //  Eigen::MatrixXd betas = as<Eigen::MatrixXd>(as<NumericMatrix>(beta_return["betas"]));
  //  Eigen::MatrixXd knots = as<Eigen::MatrixXd>(as<NumericMatrix>(beta_return["knots"]));

  for (auto j = 0; j < Y.cols(); j++) {
    if (j == 0) {
      eta.col(j) = Y.col(j).array() - betas(0, 0);
    } else {
      Eigen::MatrixXd tX = truncated_linear_cpp3(Y.col(j - 1), knots.col(j));
      eta.col(j) = Y.col(j) - tX * betas.col(j);
    }
  }

  colm = eta.colwise().mean();

  Eigen::MatrixXd centt = eta.rowwise() - eta.colwise().mean();

  tempp = (centt.adjoint() * centt) / double(centt.rows() - 1.0);

  //  return List::create(Named("col_means") = colm,
  //                      Named("cov") = tempp);
  return;
}

/// @brief Compute the Latent Slice Sampler
/// @param Y - A (n x 1) vector of observed counts
/// @param theta_cur - A (nparams x 1) vector of starting values for the MCMC
/// @param cov_eta - A list (output of compute_cov_eta)
/// @param trans_func - A list (output of compute_transform_f_lag1)
/// @param X - A (n X 1) matrix of doses
/// @param nsamples - An integer of the number of MCMC samples to run
/// @param LAM - A number loosely specifying proposal variance
/// @return A (nsamples x nparams) matrix of MCMC samples
Eigen::MatrixXd transformed_slice_sampler_cpp3(
    const Eigen::MatrixXd& Y, Eigen::VectorXd theta_cur,
    const std::function<double(Eigen::VectorXd, const Eigen::MatrixXd&, const Eigen::MatrixXd&, const Eigen::MatrixXd&, std::function<double(Eigen::VectorXd, const Eigen::MatrixXd&, const Eigen::MatrixXd&, const ptr2&)>, std::function<double(Eigen::VectorXd, const Eigen::MatrixXd&)>, const ptr2&)>
        targett,
    const Eigen::MatrixXd& priorr,
    // List cov_eta,
    Eigen::VectorXd CM, const Eigen::MatrixXd cov,
    // List trans_func,
    Eigen::MatrixXd& betas, Eigen::MatrixXd& knots, const Eigen::MatrixXd& X, int nsamples,
    double LAM,
    std::function<
        double(Eigen::VectorXd, const Eigen::MatrixXd&, const Eigen::MatrixXd&, const ptr2&)>
        lll,
    std::function<double(Eigen::VectorXd, const Eigen::MatrixXd&)> priii, const ptr2& nonlinn
) {
  // std::function<double(Eigen::VectorXd &, const Eigen::MatrixXd& ,const Eigen::VectorXd &)>
  // target= &posterior_d_cpp;

  Eigen::MatrixXd returnSamples = Eigen::MatrixXd::Zero(nsamples, theta_cur.size());

  Eigen::VectorXd Wcur(theta_cur.size());
  //  Eigen::VectorXd CM = as<Eigen::VectorXd>(as<NumericVector>(cov_eta["col_means"]));
  //  Eigen::MatrixXd cov = as<Eigen::MatrixXd>(as<NumericMatrix>(cov_eta["cov"]));
  //  // Eigen::LLT<Eigen::MatrixXd> lltOfA(cov);
  Eigen::MatrixXd Ainv = cov.llt().matrixL();
  //  // do the same for the trans_func (knots/betas)
  //  Eigen::MatrixXd betas = as<Eigen::MatrixXd>(as<NumericMatrix>(trans_func["betas"]));
  //  Eigen::MatrixXd knots = as<Eigen::MatrixXd>(as<NumericMatrix>(trans_func["knots"]));

  for (auto j = 0; j < theta_cur.size(); j++) {
    if (j == 0) {
      Wcur[j] = theta_cur[j] - betas(0, 0) - CM(0);
    } else {
      Eigen::MatrixXd tX = trunc_lin_numeric3(theta_cur[j - 1], knots.col(j));
      Wcur[j] = theta_cur(j) - (tX * betas.col(j)).value() - CM(j);
    }
  }

  // Wcur = Ainv.llt().solve(Wcur).reshaped();
  Wcur = Ainv.colPivHouseholderQr().solve(Wcur);
  Eigen::VectorXd Wsnew = Eigen::VectorXd::Zero(Wcur.size());
  // Wsnew.setZero();
  // initialize l,s,a,b, theta_new, cut, W, Wtemp
  Eigen::VectorXd l = Eigen::VectorXd::Zero(theta_cur.size());
  Eigen::VectorXd s = 0.5 * Eigen::VectorXd::Ones(theta_cur.size());
  Eigen::VectorXd a = s;
  Eigen::VectorXd b = s;
  Eigen::VectorXd theta_new = theta_cur;
  Eigen::VectorXd Wtemp = Wsnew;

  double cut = 0.0;
  double W = 0.0;
  double Wnew = 0.0;

  for (auto i = 0; i < nsamples; i++) {
    cut = targett(theta_cur, X, Y, priorr, lll, priii, nonlinn);
    cut *= -1.0;

    // W = cut - R::rexp(1.0);
    gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
    W = cut - gsl_ran_exponential(r, 1.0);
    for (auto ii = 0; ii < Wcur.size(); ii++) {
      // l[ii] = R::runif(Wcur[ii] - s[ii]/2.0, Wcur[ii]+s[ii]/2.0);
      l[ii] = gsl_ran_flat(r, Wcur[ii] - s[ii] / 2.0, Wcur[ii] + s[ii] / 2.0);

      cut = 2.0 * abs(l[ii] - Wcur[ii]);

      // s[ii] = R::rexp(LAM) + cut;
      // LAM is the rate for rexp.  gsl expects mean where mean = 1/rate.
      s[ii] = gsl_ran_exponential(r, 1 / (LAM)) + cut;
      a[ii] = l[ii] - s[ii] / 2.0;
      b[ii] = l[ii] + s[ii] / 2.0;
    }

    bool accept = true;

    while (accept) {
      for (auto ii = 0; ii < theta_cur.size(); ii++) {
        // Wsnew[ii] = R::runif(a[ii], b[ii]);
        Wsnew[ii] = gsl_ran_flat(r, a[ii], b[ii]);
      }

      Wtemp = Ainv * Wsnew;

      for (auto kk = 0; kk < theta_new.size(); kk++) {
        if (kk == 0) {
          theta_new[kk] = Wtemp[kk] + betas(0, 0) + CM[0];
        } else {
          Eigen::MatrixXd tX = trunc_lin_numeric3(theta_new[kk - 1], knots.col(kk));
          theta_new[kk] = Wtemp[kk] + (tX * betas.col(kk)).value() + CM[kk];
        }
      }

      Wnew = targett(theta_new, X, Y, priorr, lll, priii, nonlinn);
      Wnew *= -1.0;

      if (Wnew >= W) {
        accept = false;
        theta_cur = theta_new;
        Wcur = Wsnew;
      } else {
        for (auto zz = 0; zz < theta_cur.size(); zz++) {
          if (Wsnew[zz] < Wcur[zz]) {
            a[zz] = Wsnew[zz];
          } else {
            b[zz] = Wsnew[zz];
          }
        }
      }
    }

    returnSamples.row(i) = theta_cur;
  }

  return returnSamples;
}

Eigen::MatrixXd initial_slice_sampler_cpp3(
    const Eigen::MatrixXd& Y, Eigen::VectorXd Ycur,
    const std::function<double(Eigen::VectorXd, const Eigen::MatrixXd&, const Eigen::MatrixXd&, const Eigen::MatrixXd&, std::function<double(Eigen::VectorXd, const Eigen::MatrixXd&, const Eigen::MatrixXd&, const ptr2&)>, std::function<double(Eigen::VectorXd, const Eigen::MatrixXd&)>, const ptr2&)>
        targett,
    const Eigen::MatrixXd& priorr, const Eigen::MatrixXd cov, const Eigen::MatrixXd& X,
    int nsamples, double LAM,
    std::function<
        double(Eigen::VectorXd, const Eigen::MatrixXd&, const Eigen::MatrixXd&, const ptr2&)>
        lll,
    std::function<double(Eigen::VectorXd, const Eigen::MatrixXd&)> priii, const ptr2& nonlinn
) {
  Eigen::MatrixXd returnSamples = Eigen::MatrixXd::Zero(nsamples, Ycur.size());

  // Eigen::LLT<Eigen::MatrixXd> lltOfA(cov);
  Eigen::MatrixXd Ainv = cov.llt().matrixL();
  // do the same for the trans_func (knots/betas)

  // Wcur = Ainv.llt().solve(Wcur).reshaped();
  // Wcur = Ainv.colPivHouseholderQr().solve(Wcur);
  Eigen::VectorXd Ysnew = Eigen::VectorXd::Zero(Ycur.size());
  // Wsnew.setZero();
  // initialize l,s,a,b, theta_new, cut, W, Wtemp
  Eigen::VectorXd l = Eigen::VectorXd::Zero(Ycur.size());
  Eigen::VectorXd s = 0.5 * Eigen::VectorXd::Ones(Ycur.size());
  Eigen::VectorXd a = s;
  Eigen::VectorXd b = s;
  Eigen::VectorXd Ynew = Ycur;
  Eigen::VectorXd Yold = Ycur;

  double cut = 0.0;
  double W = 0.0;
  double Wnew = 0.0;

  for (auto i = 0; i < nsamples; i++) {
    // cut = as<double>(target(Ycur, X, Y));
    cut = targett(Ycur, X, Y, priorr, lll, priii, nonlinn);
    cut *= -1.0;

    // W = cut - R::rexp(1.0);
    gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
    W = cut - gsl_ran_exponential(r, 1.0);

    for (auto ii = 0; ii < Ycur.size(); ii++) {
      // l[ii] = R::runif(0.0 - s[ii]/2.0, 0.0 + s[ii]/2.0);
      l[ii] = gsl_ran_flat(r, 0.0 - s[ii] / 2.0, 0.0 + s[ii] / 2.0);

      cut = 2.0 * abs(l[ii] - 0.0);

      // s[ii] = R::rexp(LAM) + cut;
      s[ii] = gsl_ran_exponential(r, 1 / LAM) + cut;
      a[ii] = l[ii] - s[ii] / 2.0;
      b[ii] = l[ii] + s[ii] / 2.0;
    }

    bool accept = true;
    Yold = Ycur;

    while (accept) {
      for (auto ii = 0; ii < Ycur.size(); ii++) {
        // Ysnew[ii] = R::runif(a[ii], b[ii]);
        Ysnew[ii] = gsl_ran_flat(r, a[ii], b[ii]);
      }

      Ynew = Ainv * Ysnew + Ycur;

      //   Wnew = as<double>(target(Ynew, X, Y));
      Wnew = targett(Ynew, X, Y, priorr, lll, priii, nonlinn);
      Wnew *= -1.0;

      if (Wnew >= W) {
        accept = false;
        Ycur = Ynew;
        W = Wnew;
      } else {
        for (auto zz = 0; zz < Ycur.size(); zz++) {
          if (Ysnew[zz] < 0.0) {
            a[zz] = Ysnew[zz];
          } else {
            b[zz] = Ysnew[zz];
          }
        }
      }
    }

    returnSamples.row(i) = Ycur;
  }

  return returnSamples;
}

double neg_log_prior_cpp3(const Eigen::VectorXd& theta, const Eigen::MatrixXd& prior_spec) {
  double returnV = 0;
  double mean = 0;
  double sd = 0;

  const double recsqrt2pi_const = 2.0 / sqrt(pi_const);  // 2*1/sqrt(pi)
  const double recsqrt2_const = 1.0 / sqrt(2.0);         // 1/sqrt(2)
  // loop over the prior specification in prior_spec
  // when  it is 1 - Normal Prior
  // when  it is 2 - Log normal prior.
  for (int i = 0; i < theta.size(); i++) {
    int t = int(prior_spec(i, 0));

    switch (t) {
      case 1:
        mean = (theta[i] - prior_spec(i, 1));
        sd = prior_spec(i, 2);
        returnV += -log(sd) - 0.5 * mean * mean / (sd * sd) +
                   // double(theta.size()) * log(0.5 * M_2_SQRTPI * M_SQRT1_2);
                   double(theta.size()) * log(0.5 * recsqrt2pi_const * recsqrt2_const);
        break;
      case 2:

        if (theta[i] <= 0 || theta[i] > 20) {
          returnV += -std::numeric_limits<double>::infinity();
        } else {
          mean = (log(theta[i]) - prior_spec(i, 1));
          sd = prior_spec(i, 2);
          returnV += -log(sd) - log(theta[i]) - 0.5 * mean * mean / (sd * sd) +
                     // double(theta.size()) * log(0.5 * M_2_SQRTPI * M_SQRT1_2);
                     double(theta.size()) * log(0.5 * recsqrt2pi_const * recsqrt2_const);
        }
        break;
      case 3:
        returnV += log_beta_pdf(theta[i], prior_spec(i, 1), prior_spec(i, 2));
        break;
      case 4:  // gamma distribution for dirichlet
        returnV += log_gamma_pdf(theta[i], prior_spec(i, 1), prior_spec(i, 2));
        break;
      case 5: {  // Student-t: df, mu, sigma
        double df = prior_spec(i, 1);
        double mu = prior_spec(i, 2);
        double sigma = prior_spec(i, 3);
        double z = (theta[i] - mu) / sigma;

        if (sigma <= 0 || df <= 0) {
          returnV += -std::numeric_limits<double>::infinity();
        } else {
          double term1 = std::lgamma((df + 1.0) / 2.0) - std::lgamma(df / 2.0);
          // double term2 = -0.5 * log(df * M_PI) - log(sigma);
          double term2 = -0.5 * log(df * pi_const) - log(sigma);
          double term3 = -(df + 1.0) / 2.0 * log(1.0 + (z * z) / df);
          returnV += term1 + term2 + term3;
        }
        break;
      }

      default:
        // returnV -= log(0.5 * M_2_SQRTPI * M_SQRT1_2);
        returnV -= log(0.5 * recsqrt2pi_const * recsqrt2_const);
        break;
    }
  }

  return -1.0 * returnV;
}

double no_priorr(const Eigen::VectorXd& theta, const Eigen::MatrixXd& prior_spec) { return 0.0; }

////binomial likelihood
/// @brief
/// @param params
/// @param X
/// @param Y
/// @param nly
/// @return
double ll_binomial(
    Eigen::VectorXd params, const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y, const ptr2& nly
) {
  double ll = 0.0;
  Eigen::VectorXd mean = nly(params, X);

  for (int ii = 0; ii < mean.size(); ii++) {
    if (std::isnan(mean[ii])) {
      return std::numeric_limits<double>::infinity();
    }
  }
  ll += (Y.col(0).array() * mean.array().log() +
         (Y.col(1).array() - Y.col(0).array()) * (1 - mean.array()).log())
            .array()
            .sum();
  if (std::isnan(ll)) {
    ll = -std::numeric_limits<double>::infinity();
  }

  return -1.0 * ll;
}

// poisson likelihood
/// @brief
/// @param params
/// @param X
/// @param Y
/// @param nly
/// @return
double ll_poisson(
    Eigen::VectorXd params, const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y, const ptr2& nly
) {
  double ll = 0.0;

  Eigen::VectorXd log_mean = nly(params, X);

  ll += (Y.col(0).array() * log_mean.array() - exp(log_mean.array())).array().sum();
  return -1.0 * ll;
}

// random linear with added quadratic term
// double ll_continuous(Eigen::VectorXd params,
//                       const Eigen::MatrixXd& X,
//                       const Eigen::MatrixXd& Y,
//                       const ptr2& nly){
//     double ll = 0.0;
//     Eigen::VectorXd log_mean = nly(params, X);
//     double var = params[params.size() - 1.0];
//     Eigen::VectorXd sqerr = pow(Y.col(0).array() - log_mean.array(), 2.0);
//
//     ll += (-(0.5) * log(2.0 * M_PI * var) - (1 / (2.0 * var)) * sqerr.array()).sum();
//     return -1.0*ll;
// }

double ll_continuous(
    Eigen::VectorXd params, const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y, const ptr2& nly
) {
  double m0 = params[0];
  double m1 = params[1];
  double prec0 = params[params.size() - 2];
  double prec1 = params[params.size() - 1];

  double ll = 0.0;

  // compute transformed parameters
  double sigma0sq = 1.0 / prec0;
  double sigma1sq = 1.0 / prec1;
  double rho = log(sigma1sq / sigma0sq) / log(m1 / m0);
  double alpha = log(sigma0sq / pow(m0, rho));

  Eigen::VectorXd mu = nly(params, X);  // predicted means from model

  for (int i = 0; i < Y.rows(); i++) {
    double m = mu[i];
    double var = std::exp(alpha) * std::pow(m, rho);  // mean-dependent variance
    double resid = Y(i, 0) - m;
    // ll += -0.5 * std::log(2 * M_PI * var) - 0.5 * (resid * resid) / var;
    ll += -0.5 * std::log(2 * pi_const * var) - 0.5 * (resid * resid) / var;
  }

  return -1.0 * ll;
}

double ll_continuous_summary(
    const Eigen::VectorXd& params, const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y,
    const ptr2& nly
) {
  double m0 = params[0];
  double m1 = params[1];
  double prec0 = params[params.size() - 2];
  double prec1 = params[params.size() - 1];

  double sigma0sq = 1.0 / prec0;
  double sigma1sq = 1.0 / prec1;
  double rho = std::log(sigma1sq / sigma0sq) / std::log(m1 / m0);
  double alpha = std::log(sigma0sq / std::pow(m0, rho));

  Eigen::VectorXd mu = nly(params, X);  // model-predicted means

  double ll = 0.0;
  for (int i = 0; i < Y.rows(); ++i) {
    double ybar = Y(i, 0);  // observed mean
    double s = Y(i, 1);     // observed sample std
    double n = Y(i, 2);     // sample size
    double m = mu[i];       // predicted mean

    double sigma2 = std::exp(alpha) * std::pow(m, rho);
    double term1 = n * std::pow(ybar - m, 2);
    double term2 = (n - 1.0) * std::pow(s, 2);

    ll += -0.5 * n * std::log(2 * pi_const) - 0.5 * n * std::log(sigma2) -
          0.5 * (term1 + term2) / sigma2;
  }

  return -1.0 * ll;
}

// double ll_continuous_cv(Eigen::VectorXd params,
//                        const Eigen::MatrixXd& X,
//                        const Eigen::MatrixXd& Y,
//                        const ptr2& nly)
//{
//   double alpha = params[params.size() - 1];
//   double var = std::exp(alpha);
//
//   Eigen::VectorXd mu = nly(params, X);
//
//   double ll = 0.0;
//   for (int i = 0; i < Y.rows(); ++i) {
//     double resid = Y(i, 0) - mu[i];
//     ll += -0.5 * std::log(2.0 * M_PI * var) - 0.5 * (resid * resid) / var;
//   }
//
//   return -1.0 * ll;
// }

double ll_continuous_cv(
    Eigen::VectorXd params, const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y, const ptr2& nly
) {
  double prec0 = params[params.size() - 1];
  double var = 1.0 / prec0;

  Eigen::VectorXd mu = nly(params, X);

  double ll = 0.0;
  for (int i = 0; i < Y.rows(); ++i) {
    double resid = Y(i, 0) - mu[i];
    // ll += -0.5 * std::log(2.0 * M_PI * var) - 0.5 * (resid * resid) / var;
    ll += -0.5 * std::log(2.0 * pi_const * var) - 0.5 * (resid * resid) / var;
  }

  return -1.0 * ll;
}

double ll_continuous_cv_summary(
    const Eigen::VectorXd& params, const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y,
    const ptr2& nly
) {
  double prec0 = params[params.size() - 1];
  double sigma2 = 1.0 / prec0;

  Eigen::VectorXd mu = nly(params, X);  // predicted means

  double ll = 0.0;
  for (int i = 0; i < Y.rows(); ++i) {
    double ybar = Y(i, 0);
    double s = Y(i, 1);
    double n = Y(i, 2);
    double m = mu[i];

    double term1 = n * std::pow(ybar - m, 2);
    double term2 = (n - 1.0) * std::pow(s, 2);

    // ll += -0.5 * n * std::log(2 * M_PI) - 0.5 * n * std::log(sigma2) - 0.5 * (term1 + term2) /
    // sigma2;
    ll += -0.5 * n * std::log(2 * pi_const) - 0.5 * n * std::log(sigma2) -
          0.5 * (term1 + term2) / sigma2;
  }

  return -1.0 * ll;
}

double ll_lognormal_cv(
    const Eigen::VectorXd& params, const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y,
    const ptr2& nly
) {
  double prec0 = params[params.size() - 1];
  // double sdlog = std::exp(alpha);
  // double var = sdlog * sdlog;
  double var = 1.0 / prec0;

  Eigen::VectorXd mu = nly(params, X);  // this is the mean of log(Y)
  if ((mu.array().isNaN()).any()) {
    return std::numeric_limits<double>::quiet_NaN();  // or -INFINITY
  }

  double ll = 0.0;
  for (int i = 0; i < Y.rows(); ++i) {
    double logy = std::log(Y(i, 0));
    double resid = logy - mu[i];
    //    ll += -std::log(Y(i, 0))                  // Jacobian term
    //          - 0.5 * std::log(2.0 * M_PI * var)  // normalization
    //          - 0.5 * (resid * resid) / var;      // squared error
    ll += -std::log(Y(i, 0))                      // Jacobian term
          - 0.5 * std::log(2.0 * pi_const * var)  // normalization
          - 0.5 * (resid * resid) / var;          // squared error
  }

  return -1.0 * ll;  // log-likelihood; negate if needed by your sampler
}

double ll_lognormal_cv_summary(
    const Eigen::VectorXd& params, const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y,
    const ptr2& nly
) {
  double prec0 = params[params.size() - 1];
  double gamma2 = 1.0 / prec0;

  Eigen::VectorXd mu_log = nly(params, X);  // predicted log means (log[f(d_i | θ)])

  double ll = 0.0;
  for (int i = 0; i < Y.rows(); ++i) {
    double ybar = Y(i, 0);  // observed mean on original scale
    double s = Y(i, 1);     // std dev on original scale
    double n = Y(i, 2);

    double ybar_log = std::log(ybar) - 0.5 * std::log((s * s) / (ybar * ybar) + 1);
    double s_log2 = std::log((s * s) / (ybar * ybar) + 1);
    double jacobian = n * ybar_log;

    double term1 = n * std::pow(ybar_log - mu_log[i], 2);
    double term2 = (n - 1.0) * s_log2;

    //    ll += -jacobian - 0.5 * n * std::log(2 * M_PI) - 0.5 * n * std::log(gamma2) -
    //          0.5 * (term1 + term2) / gamma2;
    ll += -jacobian - 0.5 * n * std::log(2 * pi_const) - 0.5 * n * std::log(gamma2) -
          0.5 * (term1 + term2) / gamma2;
  }

  return -1.0 * ll;
}

double ll_nested_cv(
    Eigen::VectorXd params, const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y, const ptr2& nly
) {
  double m0 = params[0];
  double m1 = params[1];
  double prec0 = params[params.size() - 3];
  double prec1 = params[params.size() - 2];
  double precF = params[params.size() - 1];

  double sigma0sq = 1.0 / prec0;
  double sigma1sq = 1.0 / prec1;
  double rho = log(sigma1sq / sigma0sq) / log(m1 / m0);
  double alpha = log(sigma0sq / pow(m0, rho));

  double v_within = 1.0 / precF;

  Eigen::VectorXd mu = nly(params, X);

  double ll = 0.0;
  const int n = (int)Y.rows();
  for (int i = 0; i < n; ++i) {
    double lmean = Y(i, 0);
    double Ni = Y(i, 2);
    double s2hat = Y(i, 1);

    double m = mu[i];

    double var_between = std::exp(alpha) * std::pow(m, rho);
    double var_mean = v_within / Ni + var_between;

    double resid = lmean - m;
    // ll += -0.5 * (std::log(2.0 * M_PI) + std::log(var_mean) + (resid * resid) / var_mean);
    ll += -0.5 * (std::log(2.0 * pi_const) + std::log(var_mean) + (resid * resid) / var_mean);
  }

  return -1.0 * ll;  // negative log-likelihood
}

double ll_nested_ncv(
    Eigen::VectorXd params, const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y, const ptr2& nly
) {
  double m0 = params[0];
  double m1 = params[1];
  double prec0 = params[params.size() - 4];
  double prec1 = params[params.size() - 3];
  double precF0 = params[params.size() - 2];
  double precF1 = params[params.size() - 1];

  double sigma0sq = 1.0 / prec0;
  double sigma1sq = 1.0 / prec1;
  double sigmaF0sq = 1.0 / precF0;
  double sigmaF1sq = 1.0 / precF1;
  double rho = log(sigma1sq / sigma0sq) / log(m1 / m0);
  double alpha = log(sigma0sq / pow(m0, rho));
  double bf = log(sigmaF1sq / sigmaF0sq) / log(m1 / m0);
  double af = log(sigmaF0sq / pow(m0, bf));

  Eigen::VectorXd mu = nly(params, X);

  double ll = 0.0;
  const int n = (int)Y.rows();
  for (int i = 0; i < n; ++i) {
    double lmean = Y(i, 0);
    double Ni = Y(i, 2);
    double s2hat = Y(i, 1);

    double m = mu[i];

    double var_between = std::exp(alpha) * std::pow(m, rho);
    double v_within = std::exp(af) * std::pow(m, bf);
    double var_mean = v_within / Ni + var_between;

    double resid = lmean - m;
    // ll += -0.5 * (std::log(2.0 * M_PI) + std::log(var_mean) + (resid * resid) / var_mean);
    ll += -0.5 * (std::log(2.0 * pi_const) + std::log(var_mean) + (resid * resid) / var_mean);
  }

  return -1.0 * ll;  // negative log-likelihood
}

// generalized wrapper
double full_postt(
    Eigen::VectorXd params, const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y,
    const Eigen::MatrixXd& priorr,
    std::function<
        double(Eigen::VectorXd, const Eigen::MatrixXd&, const Eigen::MatrixXd&, const ptr2&)>
        lll,
    std::function<double(Eigen::VectorXd, const Eigen::MatrixXd&)> priorrr, const ptr2& nonlinn
) {
  double retV = lll(params, X, Y, nonlinn) + priorrr(params, priorr);

  return retV;
}

/// @brief Run the latent slice sampler
/// @param X Covariates matrix (first column is doses)
/// @param Y Counts vector
/// @param initial_val Initialization vector for the MCMC (starting point)
/// @param cov Covariance matrix for proposal
/// @param priorr Matrix of prior distributions (follows ToxicR formatting but ignoring bounds for
/// now)
/// @param model_typ Integer corresponding to model type
/// @param burnin_samples Number of samples per warm-up round
/// @param keep_samples Number of MCMC samples to keep
/// @param nrounds Number of warm-up rounds to explore for MCMC (after an initialization run)
/// @param qtiles Quantiles for the splines
/// @param LAM Changeable parameter related to proposal step size
/// @param pri_typ Integer corresponding to prior likelihood function
/// @param ll_typ Integer corresponding to the log likelihood function type
/// @return A matrix (keep_samples * # params) of MCMC samples (after burn-in)
// [[Rcpp::export]]
Eigen::MatrixXd run_latentslice_functional_general(
    const Eigen::MatrixXd& X, const Eigen::MatrixXd& Y, const Eigen::VectorXd& initial_val,
    const Eigen::MatrixXd& cov, const Eigen::MatrixXd& priorr, int model_typ, int burnin_samples,
    int keep_samples, int nrounds,
    // const Rcpp::NumericVector qtiles,
    const Eigen::VectorXd qtiles, double LAM, int pri_typ, int ll_typ
) {
  // get the nonlinearity type (1 = Hill, 2 = Power, 3 = Exponential, 5 = Logistic, otherwise use
  // Linear)
  const ptr2 nonlinn = choose_nonlinearity2(model_typ);

  // big full posterior function
  std::function<
      double(Eigen::VectorXd, const Eigen::MatrixXd&, const Eigen::MatrixXd&, const Eigen::MatrixXd&, std::function<double(Eigen::VectorXd, const Eigen::MatrixXd&, const Eigen::MatrixXd&, const ptr2&)>, std::function<double(Eigen::VectorXd, const Eigen::MatrixXd&)>, const ptr2&)>
      postt = full_postt;

  // qRcout << Y << std::endl;
  std::function<double(Eigen::VectorXd, const Eigen::MatrixXd&)> priii;
  if (pri_typ == 55) {
    priii = no_priorr;
  } else {
    priii = neg_log_prior_cpp3;
  }

  std::function<
      double(Eigen::VectorXd, const Eigen::MatrixXd&, const Eigen::MatrixXd&, const ptr2&)>
      logli;
  if (ll_typ == 55) {
    logli = ll_continuous;
  } else if (ll_typ == 56) {
    logli = ll_continuous_cv;
  } else if (ll_typ == 57) {
    logli = ll_lognormal_cv;
  } else if (ll_typ == 58) {
    logli = ll_continuous_summary;
  } else if (ll_typ == 59) {
    logli = ll_continuous_cv_summary;
  } else if (ll_typ == 60) {
    logli = ll_lognormal_cv_summary;
  } else if (ll_typ == 61) {
    logli = ll_nested_cv;
  } else if (ll_typ == 62) {
    logli = ll_nested_ncv;
  } else if (ll_typ == 66) {
    logli = ll_poisson;
  } else {
    logli = ll_binomial;
  }

  //  Eigen::MatrixXd init_samps(5,4);
  //  init_samps << 10.53683, 15.83323, 0.7422079, 0.9311202,
  //                10.51171, 15.85932, 0.7576000, 0.8718168,
  //                10.47343, 15.86369, 0.7519050, 0.8723461,
  //                10.45887, 15.68575, 0.7448072, 0.8143091,
  //                10.46783, 15.69122, 0.7583600, 0.7978510;

  //  // initial run
  Eigen::MatrixXd init_samps = initial_slice_sampler_cpp3(
      Y, initial_val, postt, priorr, cov, X, burnin_samples, LAM, logli, priii, nonlinn
  );
  Eigen::MatrixXd betas = Eigen::MatrixXd::Zero(qtiles.size() + 2, Y.cols());
  Eigen::MatrixXd knots = Eigen::MatrixXd::Zero(qtiles.size(), Y.cols());
  Eigen::VectorXd CM;
  Eigen::MatrixXd covtmp;
  for (auto ii = 0; ii < nrounds; ii++) {
    // List function_return = compute_transform_f_lag1_cpp3(init_samps, qtiles);
    compute_transform_f_lag1_cpp3(init_samps, qtiles, betas, knots);
    // List cov_eta = compute_cov_eta_cpp3(init_samps, function_return);
    compute_cov_eta_cpp3(init_samps, betas, knots, CM, covtmp);
    Eigen::VectorXd new_start = init_samps.row(burnin_samples - 1);
    //    init_samps = transformed_slice_sampler_cpp3(Y, new_start, postt, priorr,
    //                                               cov_eta, function_return, X, burnin_samples,
    //                                               LAM, logli, priii, nonlinn);
    init_samps = transformed_slice_sampler_cpp3(
        Y, new_start, postt, priorr, CM, cov, betas, knots, X, burnin_samples, LAM, logli, priii,
        nonlinn
    );
  }
  //  List function_return = compute_transform_f_lag1_cpp3(init_samps, qtiles);
  // betas.setZero();
  // knots.setZero();
  compute_transform_f_lag1_cpp3(init_samps, qtiles, betas, knots);
  //  List cov_eta = compute_cov_eta_cpp3(init_samps, function_return);
  // CM.resize(0);
  // cov.resize(0,0);
  compute_cov_eta_cpp3(init_samps, betas, knots, CM, covtmp);
  Eigen::VectorXd new_start = init_samps.row(burnin_samples - 1);
  //  init_samps = transformed_slice_sampler_cpp3(Y, new_start, postt, priorr,
  //                                             cov_eta, function_return, X, keep_samples, LAM,
  //                                             logli, priii, nonlinn);

  init_samps = transformed_slice_sampler_cpp3(
      Y, new_start, postt, priorr, CM, cov, betas, knots, X, keep_samples, LAM, logli, priii,
      nonlinn
  );
  return init_samps;
}

// double (*getLogLikeFunc(Eigen::VectorXd, const Eigen::MatrixXd&, const Eigen::MatrixXd&, const
// ptr2&))(int ll_type){
LogLikeFunction getLogLikeFunc(int ll_type) {
  std::function<
      double(Eigen::VectorXd, const Eigen::MatrixXd&, const Eigen::MatrixXd&, const ptr2&)>
      logli;

  if (ll_type == 55) {
    logli = ll_continuous;
  } else if (ll_type == 56) {
    logli = ll_continuous_cv;
  } else if (ll_type == 57) {
    logli = ll_lognormal_cv;
  } else if (ll_type == 58) {
    logli = ll_continuous_summary;
  } else if (ll_type == 59) {
    logli = ll_continuous_cv_summary;
  } else if (ll_type == 60) {
    logli = ll_lognormal_cv_summary;
  } else if (ll_type == 61) {
    logli = ll_nested_cv;
  } else if (ll_type == 62) {
    logli = ll_nested_ncv;
  } else if (ll_type == 66) {
    logli = ll_poisson;
  } else {
    logli = ll_binomial;
  }

  return logli;
}

void rg(int iter, Eigen::VectorXd mu, Eigen::MatrixXd sigma, Eigen::MatrixXd& sample) {
  int dim = sample.cols();

  Eigen::LLT<Eigen::MatrixXd> llt(sigma);
  Eigen::MatrixXd L = llt.matrixL();

  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<> d(0, 1);

  Eigen::VectorXd sample_row(dim);
  Eigen::VectorXd z(dim);

  for (int k = 0; k < iter; k++) {
    for (int i = 0; i < dim; ++i) {
      z(i) = d(gen);
    }
    sample.row(k) = mu + L * z;
  }
}

double dg(Eigen::MatrixXd V, Eigen::VectorXd mu, Eigen::MatrixXd sigma) {
  //  int k = V.size();
  //  Eigen::VectorXd diff = V - mu;
  //
  //  // Cholesky decomp
  //  Eigen::LLT<Eigen::MatrixXd> llt(sigma);
  //  if (llt.info() != Eigen::Success) {
  //    std::cout << "dg Cholesky decomp failed" << std::endl;
  //  }
  //
  //  Eigen::VectorXd y = llt.solve(diff);
  //  double quadratic_form = diff.dot(y);
  //
  //  // log determinant: Log[Sigma] = 2*sum(log(diag(L)))
  //  // double log_det = 2.0 * llt.matrixL().template
  //  // triangularView<Eigen::Lower>().diagonal().array().log().sum(); double log_det = 2.0 *
  //  // llt.matrixL().diagonal().array().log().sum();
  //
  //  //   auto lower_view = llt.matrixL().triangularView<Eigen::Lower>();
  //  //
  //  //   //log density
  //  //   double log_density = 0.5*k*log(2.0*M_PI) - 0.5 * log_det - 0.5*quadratic_form;
  //  //
  //  //   return log_density;
  //  //
  //  Eigen::MatrixXd m(3, 3);
  //  m << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  //
  //  Eigen::MatrixXd mat = llt.matrixL();
  //  auto lower_view2 = mat.triangularView<Eigen::Lower>();
  //  double log_det = lower_view2.matrixLLT().diagonal().array().log().sum();
  //
  //  //   auto lower_view = m.triangularView<Eigen::Lower>();

  return -1;

  //  int sampleSize = mu.size();
  //
  //  for (int i=0; i<V.rows(); i++){
  //
  //    Eigen::VectorXd loopRow = V.row(i);
  //
  //    Eigen::VectorXd diff = loopRow - mu;
  //    double det_cov = sigma.determinant();
  //    Eigen::MatrixXd inv_cov = sigma.inverse();
  //    double exponent = -0.5 * diff.transpose() * inv_cov * diff;
  //    //double normalization_factor = 1.0/sqrt(pow(2*M_PI, sampleSize) *det_cov);
  //    double normalization_factor = 1.0/sqrt(pow(2*M_PI, sampleSize) *det_cov);
  //
  //  }
  //  //normalization_factor*exp(exponent);
}

// double pdf_t_location_scale(double x, double mu, double sigma, double df){
double pdf_t_location_scale(double x, double df, double mu, double sigma) {
  double stand_x = (x - mu) / sigma;
  double stand_pdf_val = gsl_ran_tdist_pdf(stand_x, df);
  double location_scale_pdf_val = stand_pdf_val / sigma;

  return location_scale_pdf_val;
}
