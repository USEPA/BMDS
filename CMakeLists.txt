cmake_minimum_required(VERSION 3.1)

project(pybind11-download NONE)

include(ExternalProject)
ExternalProject_Add(
        pybind11
        PREFIX .
        GIT_REPOSITORY "https://github.com/pybind/pybind11.git"
        GIT_TAG "80dc998efced8ceb2be59756668a7e90e8bef917" # v2.10.4
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

IF(DEFINED ENV{PYTHON_LIBRARY_DIR})
  set(PYTHON_LIBRARY_DIR $ENV{PYTHON_LIBRARY_DIR})
ENDIF()

IF(DEFINED ENV{EIGEN_INCLUDE_DIR})
  set(EIGEN_INCLUDE_DIR $ENV{EIGEN_INCLUDE_DIR})
ENDIF()
#IF(DEFINED ENV{GSL_ROOT_DIR})
#  set(GSL_ROOT_DIR $ENV{GSL_ROOT_DIR})
#ENDIF()
IF(DEFINED ENV{NLOPT_LIB_DIR})
  set(NLOPT_LIB_DIR $ENV{NLOPT_LIB_DIR})
ENDIF()
IF(DEFINED ENV{CMAKE_C_COMPILER})
  set(CMAKE_C_COMPILER $ENV{CMAKE_C_COMPILER})
ENDIF()
IF(DEFINED ENV{CMAKE_CXX_COMPILER})
  set(CMAKE_CXX_COMPILER $ENV{CMAKE_CXX_COMPILER})
ENDIF()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

project(RBMDS)

include_directories("${CMAKE_SOURCE_DIR}/src/code_base" "${CMAKE_SOURCE_DIR}/src/include" ${EIGEN_INCLUDE_DIR})

file (GLOB SOURCE_FILES "${CMAKE_SOURCE_DIR}/src/code_base/*.cpp")
file (GLOB HEADER_FILES "src/include/*.h" "src/code_base/*.h" ${EIGEN_INCLUDE_DIR})
file (GLOB PYTHON_FILES "src/pybmdscpp/*.cpp" "src/pybmdscpp/*.h")

# Set up such that XCode organizes the files
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${SOURCE_FILES} ${HEADER_FILES} ${PYTHON_FILES} )

include(pybind11.cmake)
pybind11_add_module(RBMDS 
	${SOURCE_FILES}
	${HEADER_FILES}
	${PYTHON_FILES}
)

target_link_libraries(RBMDS PUBLIC)

install(TARGETS RBMDS
  COMPONENT python
  LIBRARY DESTINATION "${PYTHON_LIBRARY_DIR}"
  )
