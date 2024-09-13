test_dichotomous <- function() {
    D = c(0.0, 25.0, 75.0, 125.0, 200.0)
    Y = c(0.0, 1.0, 7.0, 15.0, 19.0)
    N = c(20.0, 20.0, 20.0, 20.0, 20.0)
    prior = c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -18.0, 0.0, 18.0, 100.0)
    result = run_bmds_dichotomous_analysis(D,Y,N,BMD_type=1, BMR=0.1, alpha=0.05, parms=2, model=3, n=5, prior=prior, prior_cols=5, degree=0)

    return(result)
}