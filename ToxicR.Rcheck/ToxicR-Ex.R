pkgname <- "ToxicR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "ToxicR-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('ToxicR')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("MAdensity_plot")
### * MAdensity_plot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: MAdensity_plot
### Title: MAdensity_plot - Create a density plot from a model averaged
###   model.
### Aliases: MAdensity_plot

### ** Examples


...
model <- ma_continuous_fit(doses,y,model_list=model_list,
                        fit_type = "mcmc",BMD_TYPE = 'sd',BMR = 1)
MAdensity_plot(model)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("MAdensity_plot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("NTP")
### * NTP

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: NTP
### Title: Williams Trend test for
### Aliases: NTP
### Keywords: datasets

### ** Examples

add(1, 1)
add(10, 1)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("NTP", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("cleveland_plot")
### * cleveland_plot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cleveland_plot
### Title: cleveland_plot - Create a Cleveland plot from a model averaged
###   model.
### Aliases: cleveland_plot

### ** Examples


...
model = ma_dichotomous_fit(D,Y,N)
cleveland_plot(model)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cleveland_plot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ma_continuous_fit")
### * ma_continuous_fit

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ma_continuous_fit
### Title: ma_continuous_fit - Fit a model averaged continuous BMD model.
### Aliases: ma_continuous_fit

### ** Examples

Hill.p <- rbind(c(481,-250.3,70,3.3),
                c(481,-250.3,40,1.3),
                c(481,-250.2,15,1.1),
                c(481,-250.3,50,4) ,
                c(10.58,9.7,70,3.5),
                c(10.58,9.7,25,3),
                c(10.58,9.7,15,2),
                c(10.58,9.7,50,4))
hill <- data.frame(a=Hill.p[,1],b=Hill.p[,2],c=Hill.p[,3],d=Hill.p[,4])
doses <- rep(c(0,6.25,12.5,25,50,100),each=10)
dosesq <- rep(c(0,6.25,12.5,25,50,100),each=30)
mean <- cont_hill_f(as.numeric(hill[2,]),doses)
y <- rinvgauss(length(mean),mean,18528.14)
model <- ma_continuous_fit(doses, y, fit_type = "laplace", BMD_TYPE = 'sd', BMR = 1)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ma_continuous_fit", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ma_dichotomous_fit")
### * ma_dichotomous_fit

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ma_dichotomous_fit
### Title: ma_dichotomous_fit - Fit a model averaged dichotomous BMD model.
### Aliases: ma_dichotomous_fit

### ** Examples


mData <- matrix(c(0, 2,50,
                  1, 2,50,
                  3, 10, 50,
                  16, 18,50,
                  32, 18,50,
                  33, 17,50),nrow=6,ncol=3,byrow=T)
D <- mData[,1]
Y <- mData[,2]
N <- mData[,3]
model = ma_dichotomous_fit(D,Y,N)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ma_dichotomous_fit", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("single_continuous_fit")
### * single_continuous_fit

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: single_continuous_fit
### Title: single_continuous_fit - Fit a single continuous BMD model.
### Aliases: single_continuous_fit

### ** Examples

M2           <- matrix(0,nrow=5,ncol=4)
colnames(M2) <- c("Dose","Resp","N","StDev")
M2[,1] <- c(0,25,50,100,200)
M2[,2] <- c(6,5.2,2.4,1.1,0.75)
M2[,3] <- c(20,20,19,20,20)
M2[,4] <- c(1.2,1.1,0.81,0.74,0.66)
model = single_continuous_fit(M2[,1,drop=F], M2[,2:4], BMD_TYPE="sd", BMR=1, ewald = T,
                             distribution = "normal",fit_type="laplace",model_type = "hill")




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("single_continuous_fit", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("single_dichotomous_fit")
### * single_dichotomous_fit

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: single_dichotomous_fit
### Title: Fit a single dichotomous dose-response model to data.
### Aliases: single_dichotomous_fit

### ** Examples

mData <- matrix(c(0, 2,50,
                  1, 2,50,
                  3, 10, 50,
                  16, 18,50,
                  32, 18,50,
                  33, 17,50),nrow=6,ncol=3,byrow=T)
D <- mData[,1]
Y <- mData[,2]
N <- mData[,3]
model = single_dichotomous_fit(D, Y, N, model_type = "hill", fit_type = "laplace")




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("single_dichotomous_fit", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
