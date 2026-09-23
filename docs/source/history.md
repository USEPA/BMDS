# History

## bmds-ui

### Version 26.1

*released on 2026-xx-xx.*

* Added LOUD Bayesian model averaging for dichotomous data
* Added LOUD Bayesian model averaging for continuous data
* Added additional continuous models for LOUD Bayesian model averaging analyses
* Added goodness of fit statistic for LOUD Bayesian model averaging
* Added download LOUD inference data object to Actions button
* Added Cochran Armitage trend test for dichotomous and cancer/multitumor data
* Added Jonckheere-Terpstra trend test for continuous data
* Added additional plotting functionality for nested dichotomous data


### Version 25.1

*Released on 2025-04-25.*

* Added Nested Dichotomous NCTR Model
* Added Rao-Scott transformation for summary Nested Dichotomous data
* Added cloning analysis action
* Only show dataset type selector for continuous data
* Enable schema migration for future backwards compatibility
* Improve database stability for BMDS Desktop
* Package versions, security updates, etc.

### Version 24.1

*Released on 2024-11-25.*

* Added BMDS Desktop mode for execution
* Added support for Multitumor modeling
* Added support for Nested Dichotomous modeling
* Added poly-k adjustment for dichotomous datasets
* Minimum `pybmds` version increased to version `24.1`
* Initial release of `bmds-ui`; project was originally forked from [shapiromatron/bmds-server](https://github.com/shapiromatron/bmds-server)

## pybmds

### Version 26.1

*released on 2026-xx-xx.*

* Added LOUD Bayesian model averaging for dichotomous data
* Added LOUD Bayesian model averaging for continuous data
* Added additional continuous models for LOUD Bayesian model averaging analyses
* Added goodness of fit statistic for LOUD Bayesian model averaging 
* Fixed p-value bug for permutation and exact Jonckheere-Terpstra trend test

### Version 25.2 

*Released on 2026-02-04.*

* Added Cochran Armitage trend test for dichotomous data
* Added Jonckheere-Terpstra trend test for continuous data
* Added additional plotting functionality for nested dichotomous data
* Changed restriction for rho parameter in non-constant variance model to allow negative values
* Added ability to count all parameters in a model automatically for the purpose of calculating AIC and p-values.
* Fixed bug in continuous hybrid BMR type
* Added validation check for instances of zero standard deviation for continuous data
* Fixed bug in Rao-Scott transformation when incidence = 0
* Changed default nested dichotomous BMR to 0.05 Extra Risk

### Version 25.1

*Released on 2025-04-25.*

* Add Nested Dichotomous NCTR Model
* Add Rao-Scott Transformation for summary Nested Dichotomous data
* Add CDF to Word Report for dichotomous bayesian model average
* Add warnings for invalid parameter settings (min, max, initial)
* Improve plotting ranges for BMDs extrapolated beyond the dose range
* Fix bug in Quantal Linear plotting
* Fix bug in Scaled Residual calculations for continuous lognormal distributions
* Remove reporting of burn-in and # of samples for bayesian modeling since pybmds uses the Laplace Approximation

### Version 24.1

*Released on 2024-11-25.*

* Added Multitumor and Multistage Cancer model
* Add Nested Dichotomous models
* Updated `bmdscore` library; integration now uses [pybind11](https://pybind11.readthedocs.io/) for direct integration instead of using a shared library. This should improve performance and stability and make it easier to develop moving forward.
* Build distributions for windows, mac, and linux using [cibuildwheel](https://cibuildwheel.pypa.io/)
* Initial release of `pybmds`; the python components of the package were originally taken from [shapiromatron/bmds](https://github.com/shapiromatron/bmds)
