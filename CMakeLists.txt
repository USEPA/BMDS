cmake_minimum_required(VERSION 3.15)

# required to find dependencies
if(NOT DEFINED ENV{VCPKG_HOST_TRIPLET})
    message(FATAL_ERROR "Environment variable 'VCPKG_HOST_TRIPLET' is not set")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(bmdscore)

set(VCPKG_HOST_TRIPLET $ENV{VCPKG_HOST_TRIPLET})

set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)
message(STATUS "pybind11 found: ${pybind11_FOUND}")

message(STATUS "vcpkg Triplet: ${VCPKG_HOST_TRIPLET}")
set(EIGEN_PATH "${CMAKE_SOURCE_DIR}/vcpkg_installed/${VCPKG_HOST_TRIPLET}/include")
message(STATUS "Eigen Path: ${EIGEN_PATH}")

set(GSL_DIR "${CMAKE_SOURCE_DIR}/vcpkg_installed/${VCPKG_HOST_TRIPLET}")
message(STATUS "GSL Path: ${GSL_DIR}")
set(GSL_INCLUDE_DIR "${GSL_DIR}/include")
set(GSL_LIBRARY_DIR "${GSL_DIR}/lib")
find_library(GSL_LIB gsl PATHS "${GSL_LIBRARY_DIR}")

set(NLOPT_DIR "${CMAKE_SOURCE_DIR}/vcpkg_installed/${VCPKG_HOST_TRIPLET}")
message(STATUS "NLOPT Path: ${NLOPT_DIR}")
set(NLOPT_DIR_INCLUDE_DIR "${NLOPT_DIR}/include")
set(NLOPT_DIR_LIBRARY_DIR "${NLOPT_DIR}/lib")
find_library(NLOPT_LIB nlopt PATHS "${NLOPT_DIR_LIBRARY_DIR}")

file (GLOB SOURCE_FILES "${CMAKE_SOURCE_DIR}/src/bmdscore/*.cpp")
file (GLOB HEADER_FILES "src/include/*.h" "src/bmdscore/*.h" ${EIGEN_DIR})
file (GLOB PYTHON_FILES "src/pybmdscpp/*.cpp" "src/pybmdscpp/*.h")
pybind11_add_module(bmdscore ${SOURCE_FILES} ${HEADER_FILES} ${PYTHON_FILES})

target_include_directories(bmdscore PRIVATE ${CMAKE_SOURCE_DIR}/src/bmdscore)
target_include_directories(bmdscore PRIVATE ${CMAKE_SOURCE_DIR}/src/include)
target_include_directories(bmdscore PRIVATE ${EIGEN_PATH})
target_include_directories(bmdscore PRIVATE ${NLOPT_DIR_INCLUDE_DIR})
target_include_directories(bmdscore PRIVATE ${GSL_INCLUDE_DIR})

target_link_libraries(bmdscore PRIVATE ${NLOPT_LIB})
target_link_libraries(bmdscore PRIVATE ${GSL_LIB})
