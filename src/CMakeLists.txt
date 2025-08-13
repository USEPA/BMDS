# This builds bmdscore + tests
cmake_minimum_required(VERSION 3.15)

# required to find dependencies
if(NOT DEFINED ENV{VCPKG_HOST_TRIPLET})
    message(FATAL_ERROR "Environment variable 'VCPKG_HOST_TRIPLET' is not set")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
IF(WIN32)
  add_compile_definitions(RBMDS_EXPORTS)
ENDIF()

project(bmdscore)

set(VCPKG_HOST_TRIPLET $ENV{VCPKG_HOST_TRIPLET})

message(STATUS "vcpkg Triplet: ${VCPKG_HOST_TRIPLET}")
set(EIGEN_PATH "${CMAKE_SOURCE_DIR}/../vcpkg_installed/${VCPKG_HOST_TRIPLET}/include")
message(STATUS "Eigen Path: ${EIGEN_PATH}")

set(GSL_DIR "${CMAKE_SOURCE_DIR}/../vcpkg_installed/${VCPKG_HOST_TRIPLET}")
message(STATUS "GSL Path: ${GSL_DIR}")
set(GSL_INCLUDE_DIR "${GSL_DIR}/include")
set(GSL_LIBRARY_DIR "${GSL_DIR}/lib")
find_library(GSL_LIB gsl PATHS "${GSL_LIBRARY_DIR}")

set(NLOPT_DIR "${CMAKE_SOURCE_DIR}/../vcpkg_installed/${VCPKG_HOST_TRIPLET}")
message(STATUS "NLOPT Path: ${NLOPT_DIR}")
set(NLOPT_DIR_INCLUDE_DIR "${NLOPT_DIR}/include")
set(NLOPT_DIR_LIBRARY_DIR "${NLOPT_DIR}/lib")
find_library(NLOPT_LIB nlopt PATHS "${NLOPT_DIR_LIBRARY_DIR}")

# ------------------------
# --- Build bmdscore -----
# ------------------------

file (GLOB SOURCE_FILES "${CMAKE_SOURCE_DIR}/bmdscore/*.cpp")
file (GLOB HEADER_FILES "${CMAKE_SOURCE_DIR}/include/*.h" "${CMAKE_SOURCE_DIR}/bmdscore/*.h")

IF(WIN32)
	#Windows build
	add_library(bmdscore SHARED ${SOURCE_FILES} ${HEADER_FILES})
ELSE()
	# Linux or MacOS
	add_library(bmdscore STATIC ${SOURCE_FILES} ${HEADER_FILES})
ENDIF(WIN32)

# make includes public so `unit_tests` can use them
target_include_directories(bmdscore PUBLIC ${CMAKE_SOURCE_DIR}/bmdscore)
target_include_directories(bmdscore PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_include_directories(bmdscore PUBLIC ${EIGEN_PATH})
target_include_directories(bmdscore PUBLIC ${NLOPT_DIR_INCLUDE_DIR})
target_include_directories(bmdscore PUBLIC ${GSL_INCLUDE_DIR})

target_link_libraries(bmdscore PRIVATE ${NLOPT_LIB})
target_link_libraries(bmdscore PRIVATE ${GSL_LIB})

# ------------------------
# --- Build Unit Tests ---
# ------------------------

file (GLOB UNIT_TEST_SOURCE "${CMAKE_SOURCE_DIR}/unitTests/*.cpp")
add_executable(unit_tests ${UNIT_TEST_SOURCE})
target_link_libraries(unit_tests bmdscore)

# Enable code coverage analysis
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")

# Add a custom target for generating coverage report
set(PROJECT_SOURCE_DIR "bmdscore")
set(TEST_EXC "'/usr/include/*'")

# Executing tests, generate coverage data, and builds HTML report
add_custom_target(run_tests_with_coverage
	COMMAND ./unit_tests
  COMMAND "here"
  COMMAND ls "${CMAKE_SOURCE_DIR}/unitTests/*"
  COMMAND ls "${CMAKE_SOURCE_DIR}/unitTests"
  COMMAND ls "${CMAKE_SOURCE_DIR}/build/CMakeFiles/unit_tests.dir/unitTests/test_main.cpp.gcno"
  COMMAND ls "${CMAKE_SOURCE_DIR}/build/CMakeFiles/bmdscore.dir/*.gcno*"
	COMMAND ${CMAKE_COMMAND} -E env
    "LCOV_EXCL_BR_LINE=:${CMAKE_SOURCE_DIR}/unitTests/*"
    gcov -b
    ${CMAKE_SOURCE_DIR}/unitTests
    ${CMAKE_SOURCE_DIR}/build/CMakeFiles/unit_tests.dir/unitTests/test_main.cpp.gcno
    ${CMAKE_SOURCE_DIR}/build/CMakeFiles/bmdscore.dir/*.gcno
	COMMAND "here2"
  COMMAND lcov
    -d ${CMAKE_SOURCE_DIR}/build/CMakeFiles/bmdscore.dir/bmdscore/ -c
    -o ${CMAKE_SOURCE_DIR}/build/coverage.info
	COMMAND "here3"
  COMMAND lcov
    -r ${CMAKE_SOURCE_DIR}/build/coverage.info ${TEST_EXC} ${EIGEN_PATH} '*gcc*'
    -o ${CMAKE_SOURCE_DIR}/build/coverage.info
    --ignore-errors unused
	COMMAND "here4"
  COMMAND genhtml
    ${CMAKE_SOURCE_DIR}/build/coverage.info
    -o ${CMAKE_SOURCE_DIR}/build/coverage
	COMMAND "here5"
  COMMENT "Running test, generating coverage data, and generating HTML report"
)
