//File:    IDPrior.h
//Purpose: Sets the priors for the binomial model  and continuous model parameters these methods are general
//         and should be able to be used across  classes without much modification
//Creator: Matt Wheeler
//Date   : 12/18/2017
//Changes: 4/13/2018 - Changed the class IDbinomPrior to IDPrior and 
//                     added typdef for IDbinomPrior.  This allows this 
//                     basic definition to be used across types of models (continuous/dichotomous etc)
//
//
//

#pragma once
#ifndef IDPriorH
#define IDPriorH

#define _USE_MATH_DEFINES

#ifdef R_COMPILATION
    //necessary things to run in R    
    #include <RcppEigen.h>
    #include <RcppGSL.h>
#else 
    #include <Eigen/Dense>
    #include <gsl/gsl_randist.h>
#endif
#include <cmath>
#include <math.h>


// IDPrior
// Creates a class where each parameter specified is given an independent
// prior.  This is done in the matrix prior_spec. The first column of this
// matrix is either 1 or 2, which indicates if the parameter is normally (1)
// or log normally (2) distributed. The second column specifies the prior mean
// and the third column specifies the prior standard deviation. 
class IDPrior{
public:

	IDPrior(IDPrior &M){
		
		prior_spec = M.prior_spec; 
	}

	// Binomial Prior Coinstructor
	IDPrior(Eigen::MatrixXd tP):prior_spec(tP){
		// TODO: do some error checking to make sure you are passing in
		// the correct size matrix and that the bounds make sense in relation
		// to the mean. 
	};
	
	void set_prior(Eigen::MatrixXd tP){
		prior_spec = tP; 
	}
	double neg_log_prior(Eigen::MatrixXd theta);
	Eigen::MatrixXd prior_mean(); 

	Eigen::MatrixXd lowerBounds() {
		return prior_spec.col(3);
	}

	Eigen::MatrixXd upperBounds() {
		return prior_spec.col(4);
	}

//private:
	// Matrix that defines the prior specifications 
	// for all of the parameters in theta
	// First Column - Type of Prior
	// Second Column - mean
	// Third Column  - Dispersion
	// Fourth Column - Lower Bound
	// Fith Column   - Upper Bound
	Eigen::MatrixXd prior_spec;
};

//set different typdefs for the type of prior continuous or binomial
//model
typedef IDPrior IDbinomPrior ; 
typedef IDPrior IDcontinuousPrior ; 



#endif
