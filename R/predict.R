#Copyright 2021  NIEHS <matt.wheeler@nih.gov>
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

.dichotomous_predict_model <- function(fit,new_doses = NULL){

    if (is.null(new_doses)){
        test_doses = fit$data[,1]
    }else{
        test_doses = new_doses
    }

    if (fit$model=="hill"){
            f <- .dich_hill_f(fit$parameters,test_doses)
          }
          if (fit$model=="gamma"){
            f <- .dich_gamma_f(fit$parameters,test_doses)
          }
          if (fit$model == "logistic"){
            f <- .dich_logist_f(fit$parameters,test_doses)
          }
          if (fit$model=="log-logistic"){
            f <- .dich_llogist_f(fit$parameters,test_doses)
          }
          if (fit$model=="probit"){
            f <- .dich_probit_f(fit$parameters,test_doses)
          }
          if (fit$model=="log-probit"){
            f<- .dich_lprobit_f(fit$parameters,test_doses)
          }
          if (fit$model=="multistage"){
            f <- .dich_multistage_f(fit$parameters,test_doses)
          }
          if (fit$model=="qlinear"){
            f<- .dich_qlinear_f(fit$parameters,test_doses)
          }
          if (fit$model=="weibull"){
            f<- .dich_weibull_f(fit$parameters,test_doses)
          }

        returnV <- list(X = test_doses, Y = f)
        return(returnV)
}

.continuous_predict_model <- function(fit, new_doses = NULL){

        data_d = fit$data

        if (ncol(data_d) == 4 ){ #sufficient statistics
            mean <- data_d[,2,drop=F]
            se   <- data_d[,4,drop=F]/sqrt(data_d[,3,drop=F])
            doses <- data_d[,1,drop=F]
            lm_fit = lm(mean ~ doses,weights = 1/(se*se))
       }else{
            Response <- data_d[,2,drop=F]
            doses = data_d[,1,drop=F]
            lm_fit = lm(Response~doses)
       }

        if (is.null(new_doses)){
             test_doses = fit$data[,1]
        }else{
             test_doses = new_doses
        }

       if (coefficients(lm_fit)[2] < 0){
         decrease = TRUE
       }else{
         decrease = FALSE
       }

        if (fit$model=="FUNL"){
             f <- cont_FUNL_f(fit$parameters,test_doses)
        }
        if (fit$model=="hill"){
             f <- cont_hill_f(fit$parameters,test_doses)
        }
        if (fit$model=="exp-3"){
             f <- cont_exp_3_f(fit$parameters,test_doses,decrease)
        }
        if (fit$model=="exp-5"){
             f <- cont_exp_5_f(fit$parameters,test_doses)
        }
        if (fit$model=="power"){
             f <- cont_power_f(fit$parameters,test_doses)
        }
        if (fit$model=="polynomial"){
             if (length(grep(": normal-ncv", tolower(fit$full_model)))>0){
                degree = length(fit$parameters) - 2
            }else{
                degree = length(fit$parameters) - 1
            }

            f <- cont_polynomial_f(fit$parameters[1:degree],test_doses)
        }

        if (grepl("Log-Normal",fit$full_model)){
            returnV <- list(X = test_doses, Y = exp(log(as.numeric(f))+ 0.5*exp(fit$parameters[length(fit$parameters)]))) #lognormal mean 
        }else{
            returnV <- list(X = test_doses, Y = as.numeric(f))
        }
        return(returnV)
}