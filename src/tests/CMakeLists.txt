cmake_minimum_required(VERSION 3.15)

# required to find dependencies
if(NOT DEFINED ENV{VCPKG_HOST_TRIPLET})
    message(FATAL_ERROR "Environment variable 'VCPKG_HOST_TRIPLET' is not set")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(test_cpp)

set(VCPKG_HOST_TRIPLET $ENV{VCPKG_HOST_TRIPLET})

message(STATUS "vcpkg Triplet: ${VCPKG_HOST_TRIPLET}")
set(EIGEN_PATH "${CMAKE_SOURCE_DIR}/../../vcpkg_installed/${VCPKG_HOST_TRIPLET}/include")
message(STATUS "Eigen Path: ${EIGEN_PATH}")

set(GSL_DIR "${CMAKE_SOURCE_DIR}/../../vcpkg_installed/${VCPKG_HOST_TRIPLET}")
message(STATUS "GSL Path: ${GSL_DIR}")
set(GSL_INCLUDE_DIR "${GSL_DIR}/include")
set(GSL_LIBRARY_DIR "${GSL_DIR}/lib")
find_library(GSL_LIB gsl NO_DEFAULT_PATH PATHS "${GSL_LIBRARY_DIR}")

set(NLOPT_DIR "${CMAKE_SOURCE_DIR}/../../vcpkg_installed/${VCPKG_HOST_TRIPLET}")
message(STATUS "NLOPT Path: ${NLOPT_DIR}")
set(NLOPT_DIR_INCLUDE_DIR "${NLOPT_DIR}/include")
set(NLOPT_DIR_LIBRARY_DIR "${NLOPT_DIR}/lib")
find_library(NLOPT_LIB nlopt NO_DEFAULT_PATH PATHS "${NLOPT_DIR_LIBRARY_DIR}")

set(CPP_BASE ../)

include_directories(${CPP_BASE}/tests/src ${CPP_BASE}/build ${CPP_BASE}/include ${CPP_BASE}/bmdscore ${EIGEN_PATH} ${NLOPT_DIR_INCLUDE_DIR} ${GSL_INCLUDE_DIR})

##Enable code coverage analysis
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
ELSE()
	#Linux or MacOS
	find_library(BMDS_LIB bmdscore HINTS ${CPP_BASE}/build)
ENDIF(WIN32)
target_link_libraries(test_cpp PRIVATE ${BMDS_LIB} ${GSL_LIB} ${NLOPT_LIB})

IF(WIN32)
	file(COPY ${CPP_BASE}/build/Release/bmdscore.dll DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE})
ELSE()
	file(COPY ${CPP_BASE}/build/libbmdscore.a DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG})
ENDIF(WIN32)
