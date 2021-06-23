#Copyright 2020  NIEHS <matt.wheeler@nih.gov>
#   
#
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software 
#and associated documentation files (the "Software"), to deal in the Software without restriction, 
#including without limitation the rights to use, copy, modify, merge, publish, distribute, 
#sublicense, and/or sell copies of the Software, and to permit persons to whom the Software 
#is furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all copies 
#or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
#INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
#PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
#HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
#CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
#OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

create_continuous_prior <- function( prior_list,model,distribution,deg=2){
  if (class(prior_list) != "BMDmodelprior"){
    stop("Prior is not of a 'BMDmodelprior' class. A probable solution is to 
          define the prior using function `create_prior_list`.")
  }
  if (!(model %in% .continuous_models )){
    stop(cat("Model Type must be one of:",.continuous_models,"\n"))
  }
  if (!(distribution %in% .continuous_distributions )){
    stop(cat("Distribution must be one of:",.continuous_distributions,"\n"))
  }
  temp = floor(deg)
  
  if ( deg < 1){
    stop("Polynomial degree must be greater than or equal to 1.")
  }
  if ( temp != deg){
    stop("Polynomial degree must be an integer.")
  }
  
  if ("hill" == model){
    p = .check_hill(prior_list,distribution) 
  }
  
  if ("FUNL" == model){
    
  }
  
  if ("exp-5" == model){
    p = .check_exp5(prior_list,distribution) 
  }
  
  if ("exp-3" == model){
    
  }
  
  if ("polynomial" == model){
    
  }
  
  if ("power" == model){
    
  }
  return(p)
}

.check_hill <- function(prior,distribution){
  #check if the normal distribution is correctly specified
  if (distribution == "normal"){
    temp = prior[[1]]
    if (nrow(temp) != 5){
      stop("Normal Hill model prior requires 5 parameters.")
    }
    if (sum(temp[,4]>temp[,5])> 0){
      stop("One of the parameter's lower bounds is greater than the upper bound.")
    }
    prior$model = "Hill Model [normal]"
    prior$parameters <- c("a","b","c","d","log(sigma^2)")
  }
  
  if (distribution == "normal-ncv"){
    temp = prior[[1]]
    if (nrow(temp) != 6){
      stop("Normal-NCV Hill model prior requires 6 parameters.")
    }
    if (sum(temp[,4]>temp[,5])> 0){
      stop("One of the parameter's lower bounds is greater than the upper bound.")
    }
    if (temp[5,5] < 0){
      stop("The prior on \rho (parameter 5) can not have a lower bound less than zero.")
    }
    prior$model = "Hill Model [normal-ncv]"
    prior$parameters <- c("a","b","c","d","rho","log(sigma^2)")
  }
  
  if (distribution == "lognormal"){
    stop("Log-Normal/Hill specification is not allowed.")
  }
  return(prior)
}

.check_exp5 <- function(prior,distribution){
  #check if the normal distribution is correctly specified
  temp = prior[[1]]
  if (distribution == "normal"){
   
    if (nrow(temp) != 5){
      stop("Normal Exponential-5  model prior requires 5 parameters.")
    }
    if (sum(temp[,4] > temp[,5])> 0){
      stop("One of the parameter's lower bounds is greater than the upper bound.")
    }
    prior$model = "Exponential-5 Model [normal]"
    prior$parameters <- c("a","b","c","d","log(sigma^2)")
  }
  
  if (distribution == "normal-ncv"){
    if (nrow(temp) != 6){
      stop("Normal Exponential-5 Hill model prior requires 6 parameters.")
    }
    if (sum(temp[,4]>temp[,5])> 0){
      stop("One of the parameter's lower bounds is greater than the upper bound.")
    }
    if (temp[5,5] < 0){
      stop("The prior on \rho (parameter 5) can not have a lower bound less than zero.")
    }
    prior$model = "Exponential-5 [normal-ncv]"
    prior$parameters <- c("a","b","c","d","rho","log(sigma^2)")
  }
  
  if (distribution == "lognormal"){
    temp = prior[[1]]
    if (nrow(temp) != 5){
      stop("Lognormal Exponential-5  model prior requires 5 parameters.")
    }
    if (sum(temp[,4]>temp[,5])> 0){
      stop("One of the parameter's lower bounds is greater than the upper bound.")
    }
    prior$model = "Exponential-5 Model [lognormal]"
    prior$parameters <- c("a","b","c","d","log(sigma^2)")
  }
  return(prior)
}