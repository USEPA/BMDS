library(actuar)

build_single_continuous_dataset <- function(){
     M2           <- matrix(0,nrow=5,ncol=4)
     colnames(M2) <- c("Dose","Resp","N","StDev")
     M2[,1] <- c(0,25,50,100,200)
     M2[,3] <- c(20,20,19,20,20)
     M2[,2] <- c(6,5.2,2.4,1.1,0.75)
     M2[,4]<-  c(1.2,1.1,0.81,0.74,0.66)
     M2
}

validate_model <- function(model, name, parameters, bmd_estimates){
     expect_equal(name, model$full_model)
     expect_equal(parameters, model$parameters, tolerance=10e-2)
     expect_equal(setNames(bmd_estimates, c("BMD", "BMDL", "BMDU")), model$bmd, tolerance=10e-2)
}

generate_validation_code <- function(AA){
     cat("\n")
     for(i in 1:(length(AA)-3)){
          cat("validate_model(", paste0("AA$Individual_Model_", i), ", ", paste0("\"", AA[[i]]$full_model, "\""), ", ", paste(list(AA[[i]]$parameters), sep=", "), ", ", paste(list(AA[[i]]$bmd), sep=", "), ")\n")
     }
}

build_ma_dataset <- function(){
     Hill.p <- rbind(c(481,-250.3,70,3.3),
                     c(481,-250.3,40,1.3),
                     c(481,-250.2,15,1.1),
                     c(481,-250.3,50,4) ,
                     c(10.58,9.7,70,3.5),
                     c(10.58,9.7,25,3),
                     c(10.58,9.7,15,2),
                     c(10.58,9.7,50,4))
     hill <- data.frame(a=Hill.p[,1],b=Hill.p[,2],
                        c=Hill.p[,3],d=Hill.p[,4])
     
     
     doses <- rep(c(0,6.25,12.5,25,50,100),each=10)
     dosesq <- rep(c(0,6.25,12.5,25,50,100),each=30)
     
     mean <- cont_hill_f(as.numeric(hill[2,]),doses)
     y <- rinvgauss(length(mean),mean,18528.14)
     return(list(doses=doses, y=y))
}