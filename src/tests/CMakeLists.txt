cmake_minimum_required(VERSION 3.5)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(test_cpp)

set(CPP_BASE ../)

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

IF(WIN32)
  find_library(GSL_LIB gsl REQUIRED HINTS ${GSL_DIR})
ELSE()
  find_package(GSL REQUIRED)
ENDIF()
find_library(NLOPT_LIB NAMES NLopt nlopt libnlopt  REQUIRED HINTS ${NLOPT_DIR})

include_directories(${CPP_BASE}/tests/src ${CPP_BASE}/build ${CPP_BASE}/include ${CPP_BASE}/bmdscore ${EIGEN_DIR} ${NLOPT_DIR}/../include ${GSL_DIR})

link_directories(${GSL_LIB_DIR} ${NLOPT_LIB_DIR})

#Enable code coverage analysis
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/../bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/../bin)

IF(WIN32)
	find_library(BMDS_LIB bmdscore HINTS ${CPP_BASE}/build/Release)
ELSE()
	find_library(BMDS_LIB bmdscore HINTS ${CPP_BASE}/build)
ENDIF()

add_executable(test_cpp src/test_cpp.cpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")

IF(WIN32)
	#Windows build
	find_library(BMDS_LIB bmdscore HINTS ${CPP_BASE}/build/Release)
	target_link_libraries(test_cpp PRIVATE ${BMDS_LIB} ${GSL_LIB} ${NLOPT_LIB})
ELSE()
	#Linux or MacOS
	find_library(BMDS_LIB bmdscore HINTS ${CPP_BASE}/build)
	target_link_libraries(test_cpp PRIVATE ${BMDS_LIB} GSL::gsl GSL::gslcblas ${NLOPT_LIB})
ENDIF(WIN32)

IF(WIN32)
	file(COPY ${CPP_BASE}/build/Release/bmdscore.dll DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE})
ELSE()
	file(COPY ${CPP_BASE}/build/libbmdscore.a DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG})
ENDIF(WIN32)
