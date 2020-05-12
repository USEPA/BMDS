#ifdef R_COMPILATION
    //necessary things to run in R    
    #include <RcppEigen.h>
    #include <RcppGSL.h>
#else 
    #include <Eigen/Dense>
#endif

#ifndef GRADIENT_FUNCTIONS_H
#define GRADIENT_FUNCTIONS_H



void gradient(Eigen::MatrixXd v, double *g, void *data, std::function<double(Eigen::MatrixXd,void*)> math_func); 

#endif
