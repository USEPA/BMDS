// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#pragma once

// we only include RcppEigen.h which pulls Rcpp.h in for us
//#include <RcppEigen.h>
//#include <math.h>
#include <cmath>
//#include <boost/math/distributions/gamma.hpp>
#include <Eigen/Dense>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

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
//using namespace Rcpp;

//////////////////////////////////////////////////////////////////////////////
/// @brief Computes the  transformation
/// @param params - A (nparams x 1) vector of Q. Linear parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the Q. Linear transformed doses
Eigen::VectorXd binomial_cqlinear(const Eigen::VectorXd& params, const Eigen::MatrixXd& X);


/// @brief Computes the  transformation
/// @param params - A (nparams x 1) vector of multistage2 parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the multistage2 transformed doses
Eigen::VectorXd binomial_cmstage2(const Eigen::VectorXd& params, 
                                 const Eigen::MatrixXd& X);



//////////////////////////////////////////////////////////////////////////////
/// @brief Computes the  transformation
/// @param params - A (nparams x 1) vector of logprobit parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the logprobit transformed doses
Eigen::VectorXd binomial_clogprobit(const Eigen::VectorXd& params, 
                                     const Eigen::MatrixXd& X);

/// @brief
/// @param params - A (nparams x 1) vector of loglogistic parameters
/// @param X 
/// @return 
Eigen::VectorXd binomial_cloglogistic(const Eigen::VectorXd& params, 
                                     const Eigen::MatrixXd& X);

/// @brief Computes the  transformation
/// @param params - A (nparams x 1) vector of Probit parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the probit transformed doses
Eigen::VectorXd binomial_cprobit(const Eigen::VectorXd& params, const Eigen::MatrixXd& X);

/// @brief Computes the  transformation
/// @param params - A (nparams x 1) vector of D. Hill parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the D. Hill transformed doses
Eigen::VectorXd binomial_chill(const Eigen::VectorXd& params, 
                               const Eigen::MatrixXd& X);

/// @brief The Weibull Model for dichotomous data
/// @param params - A (nparams x 1) vector of Weibull parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of probability 
Eigen::VectorXd binomial_cweibull(const Eigen::VectorXd& params, 
                                 const Eigen::MatrixXd& X);

/// @brief The Logistic Model for dichotomous data
/// @param params - A (nparams x 1) vector of Logistic parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of probability 
Eigen::VectorXd binomial_clogistic(const Eigen::VectorXd& params, const Eigen::MatrixXd& X);

/// @brief The Gamma Model for dichotomous data
/// @param params - A (nparams x 1) vector of Gamma parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of probability 
Eigen::VectorXd binomial_cgamma(const Eigen::VectorXd& params, 
                                 const Eigen::MatrixXd& X);

//typedef Eigen::VectorXd (*ptr)(const Eigen::VectorXd&, const Eigen::MatrixXd&);
typedef Eigen::VectorXd (*ptr2)(const Eigen::VectorXd&, const Eigen::MatrixXd&);

/// @brief 
/// @param x 
/// @param alpha 
/// @param beta 
/// @return 
double log_beta_pdf(double x,
                   double alpha,
                   double beta);

double log_gamma_pdf(double x,
                     double alpha,
                     double beta);

/// @brief Computes the  transformation
/// @param params - A (nparams x 1) vector of multistage3 parameters (DO NOT USE)
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the multistage3 transformed doses
Eigen::VectorXd binomial_mstage3(const Eigen::VectorXd& params, 
                                 const Eigen::MatrixXd& X);

/// @brief Computes the Hill transformation (DO NOT USE)
/// @param params - A (nparams x 1) vector of Hill parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the Hill transformed doses
Eigen::VectorXd hill_cpp2(const Eigen::VectorXd& params, const Eigen::MatrixXd& X);

/// @brief Computes the Power transformation (DO NOT USE)
/// @param params - A (nparams x 1) vector of Power parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the Power transformed doses
Eigen::VectorXd power_cpp2(const Eigen::VectorXd& params, const Eigen::MatrixXd& X);

/// @brief Computes the Exponential transformation (DO NOT USE)
/// @param params - A (nparams x 1) vector of Exponential parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the Exponential transformed doses
Eigen::VectorXd exponential_cpp2(const Eigen::VectorXd& params, const Eigen::MatrixXd& X);

/// @brief Computes the Linear transformation (DO NOT USE)
/// @param params - A (nparams x 1) vector of Linear parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the Linear transformed doses
Eigen::VectorXd linear_cpp2(const Eigen::VectorXd& params, const Eigen::MatrixXd& X);

/// @brief Computes the Linear transformation 
/// @param params - A (nparams x 1) vector of Linear parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the Linear transformed doses
Eigen::VectorXd continuous_linear_transform(const Eigen::VectorXd& params, const Eigen::MatrixXd& X);

/// @brief Computes the Poly2 transformation 
/// @param params - A (nparams x 1) vector of Poly2 parameters 
/// @param X - A (n X 1) matrix of doses 
/// @return A vector (n x 1) of the Poly2 transformed doses 
Eigen::VectorXd continuous_poly2_transform(const Eigen::VectorXd& params,
                                           const Eigen::MatrixXd& X);

/// @brief Computes the Power transformation 
/// @param params - A (nparams x 1) vector of Power parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the Power transformed doses
Eigen::VectorXd continuous_power_transform(const Eigen::VectorXd& params, const Eigen::MatrixXd& X);

/// @brief Computes the Exp3 transformation 
/// @param params - A (nparams x 1) vector of Exp3 parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the Exp3 transformed doses
Eigen::VectorXd continuous_exp3_transform(const Eigen::VectorXd& params, const Eigen::MatrixXd& X);

/// @brief Computes the Exp3-log transformation 
/// @param params - A (nparams x 1) vector of Exp3-log parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the Exp3-log transformed doses
Eigen::VectorXd continuous_exp3_transform_log(const Eigen::VectorXd& params, const Eigen::MatrixXd& X);

/// @brief Computes the Exp5 transformation 
/// @param params - A (nparams x 1) vector of Exp5 parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the Exp5 transformed doses
Eigen::VectorXd continuous_exp5_transform(const Eigen::VectorXd& params, const Eigen::MatrixXd& X);

/// @brief Computes the Exp5-log transformation 
/// @param params - A (nparams x 1) vector of Exp5-log parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the Exp5-log transformed doses
Eigen::VectorXd continuous_exp5_transform_log(const Eigen::VectorXd& params, const Eigen::MatrixXd& X);

/// @brief Computes the BMDS Hill transformation 
/// @param params - A (nparams x 1) vector of BMDS Hill parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the BMDS Hill transformed doses
Eigen::VectorXd continuous_hill_transform(const Eigen::VectorXd& params, const Eigen::MatrixXd& X);

/// @brief Computes the EFSA Hill transformation 
/// @param params - A (nparams x 1) vector of EFSA Hill parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the EFSA Hill transformed doses
Eigen::VectorXd continuous_hill4_efsa_transform(const Eigen::VectorXd& params, const Eigen::MatrixXd& X);

/// @brief Computes the EFSA Hill-log transformation 
/// @param params - A (nparams x 1) vector of EFSA Hill-log parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the EFSA Hill-log transformed doses
Eigen::VectorXd continuous_hill4_efsa_transform_log(const Eigen::VectorXd& params, const Eigen::MatrixXd& X);

/// @brief Computes the EFSA Inverse Exponential transformation 
/// @param params - A (nparams x 1) vector of EFSA Inverse Exponential parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the EFSA Inverse Exponential transformed doses
Eigen::VectorXd continuous_invexp_efsa_transform(const Eigen::VectorXd& params, const Eigen::MatrixXd& X);

/// @brief Computes the EFSA Inverse Exponential-log transformation 
/// @param params - A (nparams x 1) vector of EFSA Inverse Exponential-log parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the EFSA Inverse Exponential-log transformed doses
Eigen::VectorXd continuous_invexp_efsa_transform_log(const Eigen::VectorXd& params, const Eigen::MatrixXd& X);

/// @brief Computes the EFSA Lognormal transformation 
/// @param params - A (nparams x 1) vector of EFSA Lognormal parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the EFSA Lognormal transformed doses
Eigen::VectorXd continuous_lognormal_efsa_transform(const Eigen::VectorXd& params,
                                 const Eigen::MatrixXd& X);

/// @brief Computes the EFSA Lognormal-log transformation 
/// @param params - A (nparams x 1) vector of EFSA Lognormal-log parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the EFSA Lognormal-log transformed doses
Eigen::VectorXd continuous_lognormal_efsa_transform_log(const Eigen::VectorXd& params, const Eigen::MatrixXd& X);

/// @brief Computes the EFSA Gamma transformation 
/// @param params - A (nparams x 1) vector of EFSA Gamma parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the EFSA Gamma transformed doses
Eigen::VectorXd continuous_gamma_efsa_transform(const Eigen::VectorXd& params,
                             const Eigen::MatrixXd& X);

/// @brief Computes the EFSA Gamma-log transformation 
/// @param params - A (nparams x 1) vector of EFSA Gamma-log parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the EFSA Gamma-log transformed doses
Eigen::VectorXd continuous_gamma_efsa_transform_log(const Eigen::VectorXd& params, const Eigen::MatrixXd& X);

/// @brief Computes the EFSA LMS2 transformation 
/// @param params - A (nparams x 1) vector of EFSA LMS2 parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the EFSA LMS2 transformed doses
Eigen::VectorXd continuous_lms_efsa_transform(const Eigen::VectorXd& params,
                           const Eigen::MatrixXd& X);

/// @brief Computes the EFSA LMS2-log transformation 
/// @param params - A (nparams x 1) vector of EFSA LMS2-log parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the EFSA LMS2-log transformed doses
Eigen::VectorXd continuous_lms_efsa_transform_log(const Eigen::VectorXd& params, const Eigen::MatrixXd& X);

/// @brief Computes the Logistic transformation (DO NOT USE)
/// @param params - A (nparams x 1) vector of Logistic parameters
/// @param X - A (n X 1) matrix of doses
/// @return A vector (n x 1) of the Logistic transformed doses
Eigen::VectorXd logistic_cpp2(const Eigen::VectorXd& params, const Eigen::MatrixXd& X);

/// @brief Creates the pointer to the correct nonlinearity
/// @param op - An integer; 1 = Hill, 2 = Power, 3 = Exponential, 5 = Logistic, otherwise use Linear
/// @return The function pointer
ptr2 choose_nonlinearity2(int op);

// from Dirk, https://github.com/RcppCore/Rcpp/issues/967
//Rcpp::NumericVector Quantile3(Rcpp::NumericVector x, 
//                              Rcpp::NumericVector probs);
Eigen::VectorXd Quantile3(Eigen::VectorXd x, Eigen::VectorXd probs);

Eigen::MatrixXd truncated_linear_cpp3(const Eigen::VectorXd& x,
                                      const Eigen::VectorXd& knots);

Eigen::MatrixXd trunc_lin_numeric3(const double& x,
                                  const Eigen::VectorXd knots);

//List compute_transform_f_lag1_cpp3(const Eigen::MatrixXd Y,
//                                   const Rcpp::NumericVector qtiles);
void compute_transform_f_lag1_cpp3(const Eigen::MatrixXd Y, const Eigen::VectorXd qtiles,
		Eigen::MatrixXd& beta_return, Eigen::MatrixXd& knot_return);
  

//List compute_cov_eta_cpp3(Eigen::MatrixXd Y,
//                          List beta_return);
void compute_cov_eta_cpp3(Eigen::MatrixXd Y, Eigen::MatrixXd& betas, Eigen::MatrixXd& knots,
		Eigen::VectorXd& colm, Eigen::MatrixXd& tempp);

/// @brief Compute the Latent Slice Sampler
/// @param Y - A (n x 1) vector of observed counts
/// @param theta_cur - A (nparams x 1) vector of starting values for the MCMC
/// @param cov_eta - A list (output of compute_cov_eta)
/// @param trans_func - A list (output of compute_transform_f_lag1)
/// @param X - A (n X 1) matrix of doses
/// @param nsamples - An integer of the number of MCMC samples to run
/// @param LAM - A number loosely specifying proposal variance
/// @return A (nsamples x nparams) matrix of MCMC samples
Eigen::MatrixXd transformed_slice_sampler_cpp3(const Eigen::MatrixXd& Y,
                                               Eigen::VectorXd theta_cur,
                                               const std::function<double(Eigen::VectorXd, const Eigen::MatrixXd&, const Eigen::MatrixXd&, const Eigen::MatrixXd&,
                                                                          std::function<double(Eigen::VectorXd, const Eigen::MatrixXd&, const Eigen::MatrixXd&, const ptr2&)>,
                                                                          std::function<double(Eigen::VectorXd, const Eigen::MatrixXd&)>,
                                                                          const ptr2&)> targett,
                                               const Eigen::MatrixXd& priorr,
                                               //List cov_eta,
                                               Eigen::VectorXd CM,
					       const Eigen::MatrixXd cov,
                                               //List trans_func, 
                                               Eigen::MatrixXd& betas, 
                                               Eigen::MatrixXd& knots, 
                                               const Eigen::MatrixXd& X,
                                               int nsamples,
                                               double LAM,
                                               std::function<double(Eigen::VectorXd, const Eigen::MatrixXd&, const Eigen::MatrixXd&, const ptr2&)> lll,
                                               std::function<double(Eigen::VectorXd, const Eigen::MatrixXd&)> priii,
                                              const ptr2& nonlinn);
  

Eigen::MatrixXd initial_slice_sampler_cpp3(const Eigen::MatrixXd& Y,
                                          Eigen::VectorXd Ycur,
                                          const std::function<double(Eigen::VectorXd, const Eigen::MatrixXd&, const Eigen::MatrixXd&, const Eigen::MatrixXd&,
                                                            std::function<double(Eigen::VectorXd, const Eigen::MatrixXd&, const Eigen::MatrixXd&, const ptr2&)>,
                                                            std::function<double(Eigen::VectorXd, const Eigen::MatrixXd&)>,
                                                            const ptr2&)> targett,
                                          const Eigen::MatrixXd& priorr,
                                          const Eigen::MatrixXd cov,
                                          const Eigen::MatrixXd& X,
                                          int nsamples,
                                          double LAM,
                                          std::function<double(Eigen::VectorXd, const Eigen::MatrixXd&, const Eigen::MatrixXd&, const ptr2&)> lll,
                                          std::function<double(Eigen::VectorXd, const Eigen::MatrixXd&)> priii,
                                          const ptr2& nonlinn);

double neg_log_prior_cpp3(const Eigen::VectorXd& theta, 
                          const Eigen::MatrixXd& prior_spec);

double no_priorr(const Eigen::VectorXd& theta, 
                 const Eigen::MatrixXd& prior_spec );


////binomial likelihood
/// @brief 
/// @param params 
/// @param X 
/// @param Y 
/// @param nly 
/// @return 
double ll_binomial(Eigen::VectorXd params, 
                   const Eigen::MatrixXd& X, 
                   const Eigen::MatrixXd& Y, 
                   const ptr2& nly);

//poisson likelihood
/// @brief 
/// @param params 
/// @param X 
/// @param Y 
/// @param nly 
/// @return 
double ll_poisson(Eigen::VectorXd params, 
                  const Eigen::MatrixXd& X, 
                  const Eigen::MatrixXd& Y, 
                  const ptr2& nly);

//random linear with added quadratic term
//double ll_continuous(Eigen::VectorXd params, 
//                      const Eigen::MatrixXd& X, 
//                      const Eigen::MatrixXd& Y, 
//                      const ptr2& nly);
//    double ll = 0.0;
//    Eigen::VectorXd log_mean = nly(params, X);
//    double var = params[params.size() - 1.0]; 
//    Eigen::VectorXd sqerr = pow(Y.col(0).array() - log_mean.array(), 2.0);
//
//    ll += (-(0.5) * log(2.0 * M_PI * var) - (1 / (2.0 * var)) * sqerr.array()).sum();
//    return -1.0*ll;
//}

double ll_continuous(Eigen::VectorXd params, 
                     const Eigen::MatrixXd& X, 
                     const Eigen::MatrixXd& Y, 
                     const ptr2& nly);

double ll_continuous_summary(const Eigen::VectorXd& params,
                             const Eigen::MatrixXd& X,
                             const Eigen::MatrixXd& Y,
                             const ptr2& nly);

double ll_continuous_cv(Eigen::VectorXd params, 
                        const Eigen::MatrixXd& X, 
                        const Eigen::MatrixXd& Y, 
                        const ptr2& nly);

double ll_continuous_cv_summary(const Eigen::VectorXd& params,
                                const Eigen::MatrixXd& X,
                                const Eigen::MatrixXd& Y,
                                const ptr2& nly);

double ll_lognormal_cv(const Eigen::VectorXd& params,
                       const Eigen::MatrixXd& X,
                       const Eigen::MatrixXd& Y,
                       const ptr2& nly);

double ll_lognormal_cv_summary(const Eigen::VectorXd& params,
                               const Eigen::MatrixXd& X,
                               const Eigen::MatrixXd& Y,
                               const ptr2& nly);

double ll_nested_cv(Eigen::VectorXd params,
                          const Eigen::MatrixXd& X,
                          const Eigen::MatrixXd& Y,
                          const ptr2& nly);

double ll_nested_ncv(Eigen::VectorXd params,
                    const Eigen::MatrixXd& X,
                    const Eigen::MatrixXd& Y,
                    const ptr2& nly);

//generalized wrapper
double full_postt(Eigen::VectorXd params, 
                  const Eigen::MatrixXd& X, 
                  const Eigen::MatrixXd& Y, 
                  const Eigen::MatrixXd& priorr, 
                  std::function<double(Eigen::VectorXd, const Eigen::MatrixXd&, const Eigen::MatrixXd&, const ptr2&)> lll,
                  std::function<double(Eigen::VectorXd, const Eigen::MatrixXd&)> priorrr,
                  const ptr2& nonlinn);

/// @brief Run the latent slice sampler
/// @param X Covariates matrix (first column is doses)
/// @param Y Counts vector
/// @param initial_val Initialization vector for the MCMC (starting point)
/// @param cov Covariance matrix for proposal
/// @param priorr Matrix of prior distributions (follows ToxicR formatting but ignoring bounds for now)
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
Eigen::MatrixXd run_latentslice_functional_general(const Eigen::MatrixXd &X, 
                                                   const Eigen::MatrixXd& Y, 
                                                   const Eigen::VectorXd &initial_val,
                                                   const Eigen::MatrixXd& cov, 
                                                   const Eigen::MatrixXd& priorr, 
                                                   int model_typ, 
                                                   int burnin_samples,
                                                   int keep_samples, 
                                                   int nrounds, 
                                                   //const Rcpp::NumericVector qtiles, 
                                                   const Eigen::VectorXd qtiles, 
                                                   double LAM, 
                                                   int pri_typ, 
                                                   int ll_typ);
  

