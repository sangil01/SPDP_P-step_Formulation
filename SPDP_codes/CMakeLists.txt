cmake_minimum_required(VERSION 3.16)

project(SPDP_P_step_project LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(SPDP_ENABLE_ONEMKL "Enable optional Intel oneMKL backend support" OFF)

set(GUROBI_HOME "/opt/gurobi1303/linux64" CACHE PATH "Path to the Gurobi installation")
set(ORTOOLS_HOME "/opt/or-tools" CACHE PATH "Path to the OR-Tools installation")

list(PREPEND CMAKE_PREFIX_PATH "${ORTOOLS_HOME}")
find_package(ortools CONFIG REQUIRED)

if(NOT TARGET ortools::ortools)
    message(FATAL_ERROR "OR-Tools package did not define the ortools::ortools target")
endif()

message(STATUS "Using OR-Tools package: ${ortools_DIR}")

set(GUROBI_INCLUDE_DIR "${GUROBI_HOME}/include")
set(GUROBI_CPP_DIR "${GUROBI_HOME}/src/cpp")
set(GUROBI_CORE_LIB "${GUROBI_HOME}/lib/libgurobi130.so")
set(GUROBI_LIB_DIR "${GUROBI_HOME}/lib")

if(NOT EXISTS "${GUROBI_INCLUDE_DIR}/gurobi_c++.h")
    message(FATAL_ERROR "Could not find gurobi_c++.h under ${GUROBI_INCLUDE_DIR}")
endif()

if(NOT EXISTS "${GUROBI_CPP_DIR}/Env.cpp")
    message(FATAL_ERROR "Could not find the Gurobi C++ wrapper sources under ${GUROBI_CPP_DIR}")
endif()

if(NOT EXISTS "${GUROBI_CORE_LIB}")
    message(FATAL_ERROR "Could not find the Gurobi shared library under ${GUROBI_CORE_LIB}")
endif()

file(GLOB GUROBI_CPP_SOURCES "${GUROBI_CPP_DIR}/*.cpp")

if(NOT GUROBI_CPP_SOURCES)
    message(FATAL_ERROR "Could not find Gurobi C++ wrapper sources under ${GUROBI_CPP_DIR}")
endif()

message(STATUS "Using Gurobi include dir: ${GUROBI_INCLUDE_DIR}")
message(STATUS "Using Gurobi C++ source dir: ${GUROBI_CPP_DIR}")
message(STATUS "Using Gurobi core library: ${GUROBI_CORE_LIB}")

add_library(spdp_core STATIC
    src/ReadData.cpp
    src/GenMultiGraph.cpp
    src/RouteState.cpp
    src/DualFeasibleFunctions.cpp
    src/CpSatSpdpSolver.cpp
    src/CpIncumbentAdapter.cpp
    src/TwoIndexModelCore.cpp
    src/KMinComputation.cpp
    src/PstepValidInequality.cpp
    src/TwoIndexSolver.cpp
    src/InitialIncumbent.cpp
    src/PstepFormulation.cpp
    src/PstepPricing.cpp
    src/PstepMasterCG.cpp
    src/PstepBnP.cpp
    ${GUROBI_CPP_SOURCES}
)

target_include_directories(
    spdp_core
    PUBLIC include "${GUROBI_INCLUDE_DIR}" "${GUROBI_CPP_DIR}"
    PRIVATE src
)
target_link_libraries(spdp_core PUBLIC "${GUROBI_CORE_LIB}" ortools::ortools)

add_executable(SPDP_P_step src/main.cpp)
target_link_libraries(SPDP_P_step PRIVATE spdp_core)

set(SPDP_BUILD_RPATH "${GUROBI_LIB_DIR}" "${ORTOOLS_HOME}/lib")

if(SPDP_ENABLE_ONEMKL)
    if(NOT MKL_DIR AND DEFINED ENV{MKLROOT})
        set(MKL_DIR "$ENV{MKLROOT}/lib/cmake/mkl" CACHE PATH "Path to MKLConfig.cmake")
    endif()

    find_package(MKL CONFIG REQUIRED)
    message(STATUS "Intel oneMKL enabled via MKL_DIR=${MKL_DIR}")
    target_compile_definitions(spdp_core PRIVATE SPDP_HAVE_ONEMKL=1)
    target_link_libraries(spdp_core PUBLIC MKL::MKL)

    if(DEFINED MKL_ROOT)
        list(APPEND SPDP_BUILD_RPATH "${MKL_ROOT}/lib")
        get_filename_component(SPDP_ONEAPI_ROOT "${MKL_ROOT}/../.." ABSOLUTE)
        list(APPEND SPDP_BUILD_RPATH "${SPDP_ONEAPI_ROOT}/compiler/latest/lib")
    endif()
else()
    message(STATUS "Intel oneMKL backend support disabled")
endif()

set_target_properties(SPDP_P_step PROPERTIES
    BUILD_RPATH "${SPDP_BUILD_RPATH}"
)
