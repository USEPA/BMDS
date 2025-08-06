# This builds pybmds
cmake_minimum_required(VERSION 3.20)

project(pybind11-download NONE)
set(PYBIND11_FINDPYTHON OFF)

include(ExternalProject)
ExternalProject_Add(
        pybind11
        PREFIX .
        GIT_REPOSITORY "https://github.com/pybind/pybind11.git"
        GIT_TAG "ed5057ded698e305210269dafa57574ecf964483" # v3.0.0
        SOURCE_DIR "${CMAKE_BINARY_DIR}/third-party/pybind11"
        # Override default steps with no action, we just want the clone step.
        UPDATE_COMMAND ""
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
)

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

IF(WIN32)
  add_compile_definitions(RBMDS_EXPORTS)
ENDIF()

#if(NOT CMAKE_BUILD_TYPE)
#  set(CMAKE_BUILD_TYPE Release)
#endif()

#set(CMAKE_CXX_FLAGS "-O3")
#set(CMAKE_CXX_FLAGS_RELEASE "-O3")

IF(WIN32)
  find_library(GSL_LIB gsl REQUIRED HINTS ${GSL_DIR})
ELSE()
  find_package(GSL REQUIRED)
ENDIF()
find_library(NLOPT_LIB NAMES NLopt nlopt libnlopt  REQUIRED HINTS ${NLOPT_DIR})



project(bmdscore)

include_directories("${CMAKE_SOURCE_DIR}/src/bmdscore" "${CMAKE_SOURCE_DIR}/src/include" ${EIGEN_DIR} ${NLOPT_DIR}/../include ${NLOPT_DIR} ${GSL_DIR})
link_directories(${GSL_LIB_DIR} ${NLOPT_LIB_DIR})

file (GLOB SOURCE_FILES "${CMAKE_SOURCE_DIR}/src/bmdscore/*.cpp")
file (GLOB HEADER_FILES "src/include/*.h" "src/bmdscore/*.h" ${EIGEN_DIR})
file (GLOB PYTHON_FILES "src/pybmdscpp/*.cpp" "src/pybmdscpp/*.h")

# Set up such that XCode organizes the files
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${SOURCE_FILES} ${HEADER_FILES} ${PYTHON_FILES} )

include(FetchContent)

# pybind11 cmake
set(FETCHCONTENT_QUIET FALSE)

FetchContent_Declare(
  pybind11
  URL https://github.com/pybind/pybind11/archive/refs/tags/v3.0.0.tar.gz
)

FetchContent_GetProperties(pybind11)

if(NOT pybind11_POPULATED)
  FetchContent_Populate(pybind11)
  add_subdirectory(${pybind11_SOURCE_DIR} ${pybind11_BINARY_DIR})
endif()

pybind11_add_module(bmdscore
	${SOURCE_FILES}
	${HEADER_FILES}
	${PYTHON_FILES}
)

IF(WIN32)
	target_link_libraries(bmdscore PUBLIC ${GSL_LIB} ${NLOPT_LIB})
ELSE()
	target_link_libraries(bmdscore PUBLIC GSL::gsl GSL::gslcblas ${NLOPT_LIB})
ENDIF(WIN32)
