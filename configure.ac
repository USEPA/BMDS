AC_INIT([bmds_shared_object], [1.0], [matt.wheeler@nih.gov])
AM_INIT_AUTOMAKE
AC_CONFIG_MACRO_DIRS([m4])
AC_PROG_CXX
AC_PROG_CC
AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB(gsl, [gsl_ran_gaussian], [], [
        echo "Error! You need to have libgsl installed."
        exit -1
        ])
AC_CHECK_LIB(nlopt, [nlopt_get_algorithm], [], [
        echo "Error! You need to have libnlopt installed."
        exit -1
        ])


AC_ARG_VAR(EIGEN_INCLUDE,Location of Eigen include directory e.g, /usr/include/Eigen3)

LT_PREREQ([2.2])
LT_INIT()

AC_CONFIG_FILES([code_base/Makefile
		 Makefile
		 ])

AC_OUTPUT
