cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_C_COMPILER "/usr/local/apps/gcc/9.1.0/bin/gcc")
set(CMAKE_CXX_COMPILER "/usr/local/apps/gcc/9.1.0/bin/g++")

project(test_cpp)

set(CPP_BASE ../) 

#Eigen include dir
set(EXTERNAL_LIB_DIR /home/csimmo02/BMDS/RBMDS/libs)
set(EIGEN_INCLUDE_DIR ${EXTERNAL_LIB_DIR}/eigen-3.4.0)

#GSL 
set(GSL_ROOT_DIR /usr/local/apps/gsl-2.7)
find_package(GSL REQUIRED)
#set(GSL_DIR /usr/local/apps/gsl-2.7)
#set(GSL_INCLUDE_DIR ${GSL_DIR}/include)
#set(GSL_LIB_DIR ${GSL_DIR}/lib)

#NLOPT
find_library(NLOPT_LIB NAMES NLopt nlopt libnlopt  REQUIRED HINTS ${EXTERNAL_LIB_DIR}/nlopt-2.6.2/lib64)
#find_library(NLOPT_LIB NLopt REQUIRED HINTS ${EXTERNAL_LIB_DIR}/nlopt-2.6.2/out/build/x64-Debug)
#set(NLOPT_DIR ${EXTERNAL_LIB_DIR}/nlopt-2.6.2)
#set(NLOPT_INCLUDE_DIR ${NLOPT_DIR}/include)
#set(NLOPT_LIB_DIR ${NLOPT_DIR}/lib64)

#set(NLOPT_ROOT_DIR /home/csimmo02/BMDS/RBMDS/libs/nlopt-2.6.2)
#find_library(nlopt REQUIRED /home/csimmo02/BMDS/RBMDS/libs/nlopt-2.6.2  NAMES libnlopt nlopt NLopt)

#include_directories(${CPP_BASE}/build/Debug ${CPP_BASE}/include ${EIGEN_INCLUDE_DIR} ${GSL_INCLUDE_DIR} ${NLOPT_INCLUDE_DIR})
include_directories(${CPP_BASE}/tests/src ${CPP_BASE}/build ${CPP_BASE}/include ${CPP_BASE}/code_base ${EIGEN_INCLUDE_DIR} ${NLOPT_INCLUDE_DIR} ${GSL_INCLUDE_DIR})

#lib dirs
#link_directories(${GSL_LIB_DIR} ${NLOPT_LIB_DIR})
link_directories(${GSL_LIB_DIR} ${NLOPT_LIB_DIR})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/../bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/../bin)

find_library(RBMDS_LIB RBMDS HINTS ${CPP_BASE}/build)

add_executable(test_cpp src/test_cpp.cpp)

target_link_libraries(test_cpp PRIVATE ${RBMDS_LIB} GSL::gsl GSL::gslcblas ${NLOPT_LIB})

file(COPY ${CPP_BASE}/build/libRBMDS.a DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG})
