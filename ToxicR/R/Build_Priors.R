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

.parse_prior<-function(prior){
  rV <-list()
  rV$prior <- prior$prior
  
  temp_a  <- regexpr("[[][a-zA-Z]+-*[a-zA-Z]+[]]",prior$model)
  start   <- temp_a[1] + 1
  end     <- start + attr(temp_a,"match.length") - 3
  if(temp_a == -1){
    stop("Could not find a distribution for analysis.")
  }
  rV$distribution = substr(prior$model,start,end)
  rV$model = prior$mean
  return(rV)
  
}

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
  
  p = NA
  
  if ("hill" == model){
    p = .check_hill(prior_list,distribution) 
  }
  
  if ("FUNL" == model){
    p = .check_FUNLhill(prior_list,distribution) 
  }
 
  if ("exp-5" == model){
    p = .check_exp5(prior_list,distribution) 
  }
  
  if ("exp-3" == model){
    p = .check_exp3(prior_list,distribution) 
  }
  
  if ("polynomial" == model){
    p = .check_polynomial(prior_list,distribution) 
  }
  
  if ("power" == model){
    p = .check_power(prior_list,distribution) 
  }
  class(p)<- "BMD_Bayes_continuous_model"
 
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
    stop("Log-Normal/Hill specification is not presently available.")
  }
  prior$mean = .continuous_models[1]
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
      stop("Normal Exponential-5 model prior requires 6 parameters.")
    }
    if (sum(temp[,4]>temp[,5])> 0){
      stop("One of the parameter's lower bounds is greater than the upper bound.")
    }
    if (temp[5,4] < 0){
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
  prior$mean = .continuous_models[3]
  return(prior)
}

.check_power <- function(prior,distribution){
  #check if the normal distribution is correctly specified
  if (distribution == "normal"){
    temp = prior[[1]]
    if (nrow(temp) != 4){
      stop("Normal Power model prior requires 5 parameters.")
    }
    if (sum(temp[,4]>temp[,5])> 0){
      stop("One of the parameter's lower bounds is greater than the upper bound.")
    }
    if (temp[3,4] < 0){
      stop("The power parameter d (parameter 3) can not have a lower bound less than zero.")
    }
    prior$model = "Power Model [normal]"
    prior$parameters <- c("a","b","d","log(sigma^2)")
  }
  
  if (distribution == "normal-ncv"){
    temp = prior[[1]]
    if (nrow(temp) != 5){
      stop("Normal-NCV Power model prior requires 5 parameters.")
    }
    if (sum(temp[,4]> temp[,5])> 0){
      stop("One of the parameter's lower bounds is greater than the upper bound.")
    }
    if (temp[3,4] < 0){
      stop("The power parameter d (parameter 3) can not have a lower bound less than zero.")
    }
    
    if (temp[4,4] < 0){
      stop("The prior on \rho (parameter 4) can not have a lower bound less than zero.")
    }
    
    prior$model = "Power Model [normal-ncv]"
    prior$parameters <- c("a","b","d","rho","log(sigma^2)")
  }
  
  if (distribution == "lognormal"){
    stop("Log-Normal/Power specification is not presently available.")
  }
  prior$mean = .continuous_models[4]
  return(prior)
}

.check_exp3 <- function(prior,distribution){
  #check if the normal distribution is correctly specified
  temp = prior[[1]]
  if (distribution == "normal"){
    
    if (nrow(temp) != 4){
      stop("Normal Exponential-3  model prior requires 4 parameters.")
    }
    if (sum(temp[,4] > temp[,5])> 0){
      stop("One of the parameter's lower bounds is greater than the upper bound.")
    }
    
    prior$model = "Exponential-3 Model [normal]"
    prior$parameters <- c("a","b","c","d","log(sigma^2)")
  }
  
  if (distribution == "normal-ncv"){
    if (nrow(temp) != 5){
      stop("Normal Exponential-3 model prior requires 5 parameters.")
    }
    if (sum(temp[,4]>temp[,5])> 0){
      stop("One of the parameter's lower bounds is greater than the upper bound.")
    }
    if (temp[4,4] < 0){
      stop("The prior on \rho (parameter 5) can not have a lower bound less than zero.")
    }
    prior$model = "Exponential-3 [normal-ncv]"
    prior$parameters <- c("a","b","c","d","rho","log(sigma^2)")
  }
  
  if (distribution == "lognormal"){
    temp = prior[[1]]
    if (nrow(temp) != 5){
      stop("Lognormal Exponential-3  model prior requires 4 parameters.")
    }
    if (sum(temp[,4]>temp[,5])> 0){
      stop("One of the parameter's lower bounds is greater than the upper bound.")
    }
    prior$model = "Exponential-3 Model [lognormal]"
    prior$parameters <- c("a","b","c","d","log(sigma^2)")
  }
  prior$mean = .continuous_models[2]
  temp <- prior$priors
  print(temp[1:3,])
  prior$priors = matrix(NA,nrow=nrow(temp)+1,5)
  prior$priors[1:3,] = temp[1:3,]
  prior$priors[4,]   = c(1,0,1,-100,100)
  prior$priors[5:nrow(prior$priors), ] = temp[4:nrow(temp),]
  cat("NOTE: Parameter 'c' added to prior list. It is not used in the analysis.")
  return(prior)
}

.check_polynomial <- function(prior,distribution){
  #check if the normal distribution is correctly specified
  temp = prior[[1]]
   if (sum(temp[,4]>temp[,5])> 0){
    stop("One of the parameter's lower bounds is greater than the upper bound.")
  }
  
  temp_p <- c("b0")
  for (ii in 2:(nrow(temp)-1)){
    temp_p <- c(temp_p,sprintf("b%s",ii-1))
  }
  
  if (distribution == "normal"){
    if (nrow(temp) < 3){
      stop("Normal Polynomial models require 3 or more parameters.")
    }
    prior$model = "Hill Model [normal]"
    prior$parameters <- c(temp_p,"log(sigma^2)")
  }
  
  if (distribution == "normal-ncv"){
    temp = prior[[1]]
    if (nrow(temp) < 4){
      stop("Normal Polynomial models require 4 or more parameters.")
    }

    if (temp[nrow(temp)-2,5] < 0){
      stop("The prior on \rho (parameter 5) can not have a lower bound less than zero.")
    }
    prior$model = "Hill Model [normal-ncv]"
    temp_p[length(temp_p)] = "rho"
    prior$parameters <- c(temp_p,"log(sigma^2)")
  }
  
  if (distribution == "lognormal"){
    stop("Log-Normal/Polynomial specification is not presently available.")
  }
  prior$mean = .continuous_models[6]
  return(prior)
}

.check_FUNLhill <- function(prior,distribution){
  #check if the normal distribution is correctly specified
  temp = prior[[1]]
  if (sum(temp[,4]>temp[,5])> 0){
    stop("One of the parameter's lower bounds is greater than the upper bound.")
  }
  
  if (distribution == "normal"){
    
    if (nrow(temp) != 7){
      stop("Normal FUNL model prior requires 7 parameters.")
    }
  
    prior$model = "FUNL Model [normal]"
    prior$parameters <- c("a","b","LM","LD","NM","ND","log(sigma^2)")
  }
  
  if (distribution == "normal-ncv"){

    if (nrow(temp) != 8){
      stop("Normal-NCV Hill model prior requires 8 parameters.")
    }
    
    if (temp[7,5] < 0){ #check rho
      stop("The prior on \rho (parameter 7) can not have a lower bound less than zero.")
    }
    prior$model = "FUNL Model [normal-ncv]"
    prior$parameters <- c("a","b","LM","LD","NM","ND","rho","log(sigma^2)")
  }
  
  if (distribution == "lognormal"){
    stop("Log-Normal/FUNL specification is not presently available.")
  }
  prior$mean = .continuous_models[5]
  return(prior)
}


