cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_C_COMPILER "/usr/local/apps/gcc/9.1.0/bin/gcc")
set(CMAKE_CXX_COMPILER "/usr/local/apps/gcc/9.1.0/bin/g++")
set(PYTHON_LIBRARY_DIR "/home/csimmo02/BMDS/RBMDS/RBMDS_CMAKE/venv/Lib/site-packages")

project(RBMDS VERSION 0.1.0)

#Eigen include dir
set(EXTERNAL_LIB_DIR "/home/csimmo02/BMDS/RBMDS/libs")
set(EIGEN_INCLUDE_DIR ${EXTERNAL_LIB_DIR}/eigen-3.4.0)

#GSL 
#find_library(GSL_LIB gsl REQUIRED )
#find_package(GSL REQUIRED HINTS /usr/local/apps/gsl-2.7)
#set(GSL_INCLUDE_DIR ${EXTERNAL_LIB_DIR}/gsl-master/out/build/x64-Release)

set(GSL_ROOT_DIR /usr/local/apps/gsl-2.7)
find_package(GSL REQUIRED)

#set(GSL_DIR /usr/local/apps/gsl-2.7)
#set(GSL_INCLUDE_DIR ${GSL_DIR}/include)
#set(GSL_LIB_DIR ${GSL_DIR}/lib)

#NLOPT
#find_library(NLOPT_LIB NLopt REQUIRED HINTS ${EXTERNAL_LIB_DIR}/nlopt-2.6.2)
#set(NLOPT_INCLUDE_DIR ${EXTERNAL_LIB_DIR}/nlopt-2.6.2/out/build/x64-Debug ${EXTERNAL_LIB_DIR}/nlopt-2.6.2/src/api)
#set(NLOPT_DIR ${EXTERNAL_LIB_DIR}/nlopt-2.6.2)
#set(NLOPT_INCLUDE_DIR ${NLOPT_DIR}/include)
#set(NLOPT_LIB_DIR ${NLOPT_DIR}/lib64)

find_library(NLOPT_LIB NAMES NLopt nlopt libnlopt  REQUIRED HINTS ${EXTERNAL_LIB_DIR}/nlopt-2.6.2/lib64)

# Include dir
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/build/Debug ${EIGEN_INCLUDE_DIR} ${GSL_INCLUDE_DIR} ${NLOPT_INCLUDE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/code_base ${CMAKE_CURRENT}/include ${EIGEN_INCLUDE_DIR} ${GSL_INCLUDE_DIR} $NLOPT_INCLUDE_DIR})

#lib dirs
link_directories(${GSL_LIB_DIR} ${NLOPT_LIB_DIR})

# Src
AUX_SOURCE_DIRECTORY(code_base SRC_FILES)

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
add_library(RBMDS STATIC ${SOURCE_FILES})

# Libs
#target_link_libraries(RBMDS PRIVATE ${GSL_LIB} ${NLOPT_LIB})
target_link_libraries(RBMDS PRIVATE GSL::gsl GSL::gslcblas ${NLOPT_LIB})

# Include directories
target_include_directories(RBMDS PRIVATE include)

# Install
#install(TARGETS RBMDS DESTINATION )

# Install the headers
#install(FILES include/RBMDS DESTINATION include)

# Create base directory
#install(DIRECTORY include DESTINATION include)
