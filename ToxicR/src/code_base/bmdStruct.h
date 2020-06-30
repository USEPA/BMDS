
//necessary things to run in R  
#ifndef BMD_ANALYSIS_h
#define BMD_ANALYSIS_h


enum est_method {est_mle = 1, est_laplace=2, est_mcmc=3}; 
enum dich_model {d_hill =1, d_gamma=2,d_logistic=3, d_loglogistic=4,
                 d_logprobit=5, d_multistage=6,d_probit=7,
                 d_qlinear=8,d_weibull=9}; 

enum cont_model {hill = 6,exp_3 = 3,exp_5=5,power=8, gain_loss_model = 10, polynomial = 666}; 
enum distribution {normal = 1, normal_ncv = 2, log_normal = 3}; 

// Dichotomous Structures
struct dichotomous_analysis{
  dich_model model; 
  int n; 
  double *Y; // observed +
  double *doses; // 
  double *n_group; //
  double *prior; // a column order matrix px5 where p is 
                 // the number of parametersd
  int BMD_type;  
  double BMR; 
  double alpha; 
  int degree;  // used for multistage
  int samples; // number of MCMC samples. 
  int burnin;  // size of burin
  int parms;   // number of parameters
  int prior_cols; // colunns in the prior
};

struct dichotomous_model_result{
  int      model;           // dichotomous model specification
  int      nparms; 		      //number of parameters in the model
  double  *parms;           // Parameter Estimate 
  double  *cov;             // Covariance Estimate
  double   max;             // Value of the Likelihood/Posterior at the maximum
  int      dist_numE;       // number of entries in rows for the bmd_dist
  double  *bmd_dist;        // bmd distribution (dist_numE x 2) matrix
};

struct dichotomousMA_analysis{
  int    nmodels;         //number of models for each 
  double **priors;     //priors
  int    *nparms;      //parameter in each model
  int    *actual_parms;//actual number of parameters in the model
  int    *prior_cols;  // columns in the prior if there are 'more' in the future
  int    *models;      // given model
  double *modelPriors; // prior probability on the model
};

struct dichotomousMA_result{
  int                      nmodels; //number of models for each 
  dichotomous_model_result **models; //priors
  int                    dist_numE; // number of entries in rows for the bmd_dist
  double               *post_probs; // posterior probabilities
  double                 *bmd_dist; // bmd ma distribution (dist_numE x 2) matrix
};

// Continuous Structures
struct continuous_analysis{
  cont_model model; 
  int n; 
  bool suff_stat; //true if the data are in sufficient statistics format
  double *Y; // observed +
  double *doses; // 
  double *sd; //
  double *n_group; //
  double *prior; // a column order matrix px5 where p is 
                // the number of parametersd
  int BMD_type; 
  
  bool isIncreasing; 
  double BMR; 
  double tail_prob; 
  int    disttype; 
  double alpha; 
  est_method fitter; 
  int samples; // number of MCMC samples. 
  int burnin;  // burn in 
  int parms; // number of parameters 
  int prior_cols; 
};
struct continuousMA_analysis{
  int    nmodels;         //number of models for each 
  double **priors;     //priors
  int    *nparms;      //parameter in each model
  int    *actual_parms;//actual number of parameters in the model
  int    *prior_cols;  // columns in the prior if there are 'more' in the future
  int    *models;      // given model
  int    *disttype;    // given distribution type
  double *modelPriors; // prior probability on the model
};
struct continuous_model_result{
  int      model;           // continuous model specification
  int      dist;            // distribution_type 
  int      nparms; 		      //number of parameters in the model
  double  *parms;           // Parameter Estimate 
  double  *cov;             // Covariance Estimate
  double   max;             // Value of the Likelihood/Posterior at the maximum
  int      dist_numE;       // number of entries in rows for the bmd_dist
  double  *bmd_dist;        // bmd distribution (dist_numE x 2) matrix
};
struct continuousMA_result{
  int                      nmodels; //number of models for each 
  continuous_model_result **models; //priors
  int                    dist_numE; // number of entries in rows for the bmd_dist
  double               *post_probs; // posterior probabilities
  double                 *bmd_dist; // bmd ma distribution (dist_numE x 2) matrix
};

// mcmc structures
struct bmd_analysis_MCMC{
  int model; 
  unsigned int burnin; 
  unsigned int samples;
  unsigned int nparms; 
  double * BMDS; 
  double * parms; 
};
struct ma_MCMCfits{
  unsigned int nfits; 
  bmd_analysis_MCMC **analyses;   
}; 

// odds and ends
bmd_analysis_MCMC * new_mcmc_analysis(int model,
                                       int parms, 
                                       unsigned int samples);
void del_mcmc_analysis(bmd_analysis_MCMC *an); 
void cp_prior(Eigen::MatrixXd temp, double *priors);
void del_continuousMA_analysis(continuousMA_analysis &CMA);
void del_continuous_analysis(continuous_analysis a); 
continuous_model_result * new_continuous_model_result(int model,
													  unsigned int n_parm,
                            unsigned int n_elm);
                                                      
void del_continuous_model_result(continuous_model_result * cm);   


#endif
