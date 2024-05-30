# This builds bmdscore + tests
cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

IF(DEFINED ENV{EIGEN_DIR})
  set(EIGEN_DIR $ENV{EIGEN_DIR})
ENDIF()
IF(DEFINED ENV{GSL_DIR})
  set(GSL_DIR $ENV{GSL_DIR})
  IF(NOT WIN32)
    set(GSL_ROOT_DIR $ENV{GSL_DIR})
  ENDIF()
ENDIF()
IF(DEFINED ENV{NLOPT_DIR})
  set(NLOPT_DIR $ENV{NLOPT_DIR})
ENDIF()
IF(DEFINED ENV{CMAKE_C_COMPILER})
  set(CMAKE_C_COMPILER $ENV{CMAKE_C_COMPILER})
ENDIF()
IF(DEFINED ENV{CMAKE_CXX_COMPILER})
  set(CMAKE_CXX_COMPILER $ENV{CMAKE_CXX_COMPILER})
ENDIF()
IF(DEFINED ENV{CMAKE_CXX_COMPILER})
  set(CMAKE_CXX_COMPILER $ENV{CMAKE_CXX_COMPILER})
ENDIF()
IF(DEFINED ENV{TEST_EXC})
  set(TEST_EXC $ENV{TEST_EXC})
ENDIF()

IF(WIN32)
  add_compile_definitions(RBMDS_EXPORTS)
ENDIF()

project(bmdscore VERSION 0.1.0)

#Find libraries
IF(WIN32)
  find_library(GSL_LIB gsl REQUIRED HINTS ${GSL_DIR})
ELSE()
  find_package(GSL REQUIRED)
ENDIF()
find_library(NLOPT_LIB NAMES NLopt nlopt libnlopt  REQUIRED HINTS ${NLOPT_DIR})


#lib dirs
link_directories(${GSL_LIB_DIR} ${NLOPT_LIB_DIR})

# Include dir
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/code_base ${CMAKE_CURRENT_SOURCE_DIR}/include ${EIGEN_DIR} ${GSL_DIR} ${NLOPT_DIR} ${NLOPT_DIR}/../include)

# Headers
set(PROJECT_SOURCE_DIR "code_base")
set(PROJECT_INCLUDE_DIR "include")

# Source files
set(SOURCE_FILES
        ${PROJECT_SOURCE_DIR}/analysis_of_deviance.cpp
	${PROJECT_SOURCE_DIR}/analysis_of_deviance.h
        ${PROJECT_SOURCE_DIR}/bmds_helper.cpp
	${PROJECT_SOURCE_DIR}/bmds_helper.h
	${PROJECT_SOURCE_DIR}/bmdStruct.cpp
	${PROJECT_SOURCE_DIR}/bmdStruct.h
	${PROJECT_SOURCE_DIR}/continuous_clean_aux.cpp
	${PROJECT_SOURCE_DIR}/continuous_clean_aux.h
	${PROJECT_SOURCE_DIR}/continuous_entry_code.cpp
	${PROJECT_SOURCE_DIR}/continuous_entry_code.h
	${PROJECT_SOURCE_DIR}/continuous_model_functions.cpp
	${PROJECT_SOURCE_DIR}/continuous_model_functions.h
	${PROJECT_SOURCE_DIR}/DichGammaBMD_NC.cpp
	${PROJECT_SOURCE_DIR}/DichLogisticBMD_NC.cpp
	${PROJECT_SOURCE_DIR}/DichLogProbitBMD_NC.cpp
	${PROJECT_SOURCE_DIR}/dichotomous_entry_code.cpp
	${PROJECT_SOURCE_DIR}/dichotomous_entry_code.h
	${PROJECT_SOURCE_DIR}/DichProbitBMD_NC.cpp
	${PROJECT_SOURCE_DIR}/DichQlinearBMD_NC.cpp
	${PROJECT_SOURCE_DIR}/DichWeibullBMD_NC.cpp
	${PROJECT_SOURCE_DIR}/framework.h
	${PROJECT_SOURCE_DIR}/gradient.cpp
	${PROJECT_SOURCE_DIR}/gradient.h
	${PROJECT_SOURCE_DIR}/helperfunctions.cpp
	${PROJECT_SOURCE_DIR}/IDPrior.cpp
	${PROJECT_SOURCE_DIR}/IDPriorMCMC.cpp
	${PROJECT_SOURCE_DIR}/list_r_conversion.h
	${PROJECT_SOURCE_DIR}/lognormal_EXP_NC.cpp
	${PROJECT_SOURCE_DIR}/lognormal_HILL_NC.cpp
	${PROJECT_SOURCE_DIR}/lognormal_POLYNOMIAL_NC.cpp
	${PROJECT_SOURCE_DIR}/lognormal_POWER_NC.cpp
	${PROJECT_SOURCE_DIR}/lognormalModels.cpp
        ${PROJECT_SOURCE_DIR}/normal_EXP_NC.cpp
	${PROJECT_SOURCE_DIR}/normal_FUNL_NC.cpp
        ${PROJECT_SOURCE_DIR}/normal_HILL_NC.cpp
	${PROJECT_SOURCE_DIR}/normal_likelihoods.cpp
	${PROJECT_SOURCE_DIR}/normal_POLYNOMIAL_NC.cpp
	${PROJECT_SOURCE_DIR}/normal_POWER_NC.cpp
	${PROJECT_SOURCE_DIR}/normalModels.cpp
#	${PROJECT_SOURCE_DIR}/pch.cpp
#	${PROJECT_SOURCE_DIR}/pch.h
	${PROJECT_SOURCE_DIR}/stdafx.cpp
	${PROJECT_SOURCE_DIR}/stdafx.h
	${PROJECT_SOURCE_DIR}/targetver.h
	${PROJECT_INCLUDE_DIR}/binomialTests.h
	${PROJECT_INCLUDE_DIR}/binomModels.h
	${PROJECT_INCLUDE_DIR}/bmd_calculate.h
	${PROJECT_INCLUDE_DIR}/bmds_entry.h
	${PROJECT_INCLUDE_DIR}/cBMDstatmod.h
	${PROJECT_INCLUDE_DIR}/cmodeldefs.h
	${PROJECT_INCLUDE_DIR}/dBMDstatmod.h
	${PROJECT_INCLUDE_DIR}/DichGammaBMD_NC.h
	${PROJECT_INCLUDE_DIR}/DichHillBMD_NC.h
	${PROJECT_INCLUDE_DIR}/DichLogisticBMD_NC.h
	${PROJECT_INCLUDE_DIR}/DichLogLogisticBMD_NC.h
	${PROJECT_INCLUDE_DIR}/DichLogProbitBMD_NC.h
	${PROJECT_INCLUDE_DIR}/DichMultistageBMD_NC.h
	${PROJECT_INCLUDE_DIR}/DichProbitBMD_NC.h
	${PROJECT_INCLUDE_DIR}/DichQlinearBMD_NC.h
	${PROJECT_INCLUDE_DIR}/DichWeibullBMD_NC.h
	${PROJECT_INCLUDE_DIR}/IDPrior.h
	${PROJECT_INCLUDE_DIR}/IDPriorMCMC.h
	${PROJECT_INCLUDE_DIR}/log_likelihoods.h
	${PROJECT_INCLUDE_DIR}/lognormal_EXP_NC.h
	${PROJECT_INCLUDE_DIR}/lognormal_HILL_NC.h
	${PROJECT_INCLUDE_DIR}/lognormal_likelihoods.h
	${PROJECT_INCLUDE_DIR}/lognormal_POLYNOMIAL_NC.h
	${PROJECT_INCLUDE_DIR}/lognormal_POWER_NC.h
	${PROJECT_INCLUDE_DIR}/lognormalModels.h
	${PROJECT_INCLUDE_DIR}/lognormalTests.h
	${PROJECT_INCLUDE_DIR}/mcmc_anal.h
	${PROJECT_INCLUDE_DIR}/mcmc_analysis.h
	${PROJECT_INCLUDE_DIR}/mcmc_struct.h
	${PROJECT_INCLUDE_DIR}/normal_EXP_NC.h
	${PROJECT_INCLUDE_DIR}/normal_FUNL_NC.h
	${PROJECT_INCLUDE_DIR}/normal_HILL_NC.h
	${PROJECT_INCLUDE_DIR}/normal_likelihoods.h
	${PROJECT_INCLUDE_DIR}/normal_POLYNOMIAL_NC.h
	${PROJECT_INCLUDE_DIR}/normal_POWER_NC.h
	${PROJECT_INCLUDE_DIR}/normalModels.h
	${PROJECT_INCLUDE_DIR}/statisticalmodel.h
	${PROJECT_INCLUDE_DIR}/statmod.h
	)

# Set up such that XCode organizes the files correctly
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${SOURCE_FILES})

# Add library
IF(WIN32)
	#Windows build
	add_library(bmdscore SHARED ${SOURCE_FILES})
	target_link_libraries(bmdscore PRIVATE ${GSL_LIB} ${NLOPT_LIB})
ELSE()
	#Linux or MacOS
	add_library(bmdscore STATIC ${SOURCE_FILES})
	target_link_libraries(bmdscore PRIVATE GSL::gsl GSL::gslcblas ${NLOPT_LIB})
ENDIF(WIN32)

# Include directories
target_include_directories(bmdscore PRIVATE include)

#Add the executable target for the tests
add_executable(unit_tests ./unitTests/test_main.cpp ./unitTests/unit_tests.cpp ./unitTests/integration_tests.cpp ./unitTests/test_func.cpp)

#Link the library to the test executable
target_link_libraries(unit_tests bmdscore)

#ENable code coverage analysis
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")

#Add a custom target for generating coverage report
add_custom_target(coverage
	COMMAND ${CMAKE_COMMAND} -E env "LCOV_EXCL_BR_LINE=:${PROJECT_SOURCE_DIR}/../unitTests/*" gcov -b ${PROJECT_SOURCE_DIR}../unitTests ${PROJECT_SOURCE_DIR}/../build/CMakeFiles/unit_tests.dir/unitTests/test_main.cpp.gcno ${PROJECT_SOURCE_DIR}/../build/CMakeFiles/bmdscore.dir/*.gcno

	COMMAND lcov -d ./CMakeFiles/bmdscore.dir/code_base/ -c -o ./coverage.info
	COMMAND lcov -r ./coverage.info ${TEST_EXC} ${EIGEN_DIR} '*gcc*' -o ./coverage.info --ignore-errors unused
	COMMAND lcov -r ./coverage.info
	COMMAND lcov -r ./coverage.info -o ./coverage.info
	COMMAND genhtml ./coverage.info -o ./coverage
	WORKING_DIRECTORY ./
	COMMENT "Generating code coverage report"
	message(STATUS "NLOPT_DIR=${NLOPT_DIR}")

)

#Add a custom target for executing tests and generating coverage report
add_custom_target(run_tests
	COMMAND ./unit_tests
	DEPENDS unit_tests
	WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
	COMMENT "Running tests and generating code coverage data"
)

#Add a custom target for executing tests, generating coverage data, and generating HTML report
add_custom_target(run_tests_with_coverage
	COMMAND ./unit_tests
	COMMAND ${CMAKE_COMMAND} -E env "LCOV_EXCL_BR_LINE=:${PROJECT_SOURCE_DIR}/../unitTests/*" gcov -b ${PROJECT_SOURCE_DIR}../unitTests ${PROJECT_SOURCE_DIR}/../build/CMakeFiles/unit_tests.dir/unitTests/test_main.cpp.gcno ${PROJECT_SOURCE_DIR}/../build/CMakeFiles/bmdscore.dir/*.gcno
	COMMAND lcov -c -o ./coverage.info
	COMMAND lcov -r ./coverage.info ${TEST_EXC} ${EIGEN_DIR} '/usr/include/*' '*gcc*' -o ./coverage.info --ignore-errors unused
	COMMAND genhtml ./coverage.info -o ./coverage
	WORKING_DIRECTORY ./  #${PROJECT_BINARY_DIR}
	COMMENT "Running test, generating coverage data, and generating HTML report"
)

