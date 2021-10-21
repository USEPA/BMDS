context("Single Continuous Models MCMC")

test_that("Normal Ewald Hill", {
     set.seed(5981)
     M2 <- build_single_continuous_dataset()
     c = single_continuous_fit(M2[,1,drop=F],M2[,2:4],BMD_TYPE="sd",BMR=1, ewald = T,
                               distribution = "normal",fit_type="mcmc",model_type = "hill",degree = 4)
     #FIXME MCMC has different class structure
     #expect_equal("hill", c$full_model)
     #expect_equal(c(6.07, -5.33, 40.066, 3.108, -0.189), c$parameters, tolerance=10e-2)
     expect_equal(setNames(c(24.5, 20, 29.2), c("BMD", "BMDL", "BMDU")), c$bmd, tolerance=10e-2)
})