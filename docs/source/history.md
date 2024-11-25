# History

## bmds-ui

### Version 24.1

*Released on 2024-11-25.*

* Added BMDS Desktop mode for execution
* Added support for Multitumor modeling
* Added support for Nested Dichotomous modeling
* Added poly-k adjustment for dichotomous datasets
* Minimum `pybmds` version increased to version `24.1`
* Initial release of `bmds-ui`; project was originally forked from [shapiromatron/bmds-server](https://github.com/shapiromatron/bmds-server)

## pybmds

### Version 24.1

*Released on 2024-11-25.*

* Added Multitumor and Multistage Cancer model
* Add Nested Dichotomous models
* Updated `bmdscore` library; integration now uses [pybind11](https://pybind11.readthedocs.io/) for direct integration instead of using a shared library. This should improve performance and stability and make it easier to develop moving forward.
* Build distributions for windows, mac, and linux using [cibuildwheel](https://cibuildwheel.pypa.io/)
* Initial release of `pybmds`; the python components of the package were originally taken from [shapiromatron/bmds](https://github.com/shapiromatron/bmds)
