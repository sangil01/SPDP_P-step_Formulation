cmake_minimum_required(VERSION 3.16)

project(SPDP_P_step_project LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(SPDP_ENABLE_ONEMKL "Enable optional Intel oneMKL backend support" OFF)
option(SPDP_ENABLE_LONG_CP_TESTS "Register long CP-SAT benchmark tests" OFF)

set(GUROBI_HOME "/opt/gurobi1301/linux64" CACHE PATH "Path to the Gurobi installation")
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

include(CTest)
if(BUILD_TESTING)
    add_executable(
        SPDP_cp_tests
        tests/TestMain.cpp
        tests/RouteStateTests.cpp
        tests/InitialIncumbentSearchTests.cpp
        tests/CpSatSpdpSolverTests.cpp
        tests/CpIncumbentAdapterTests.cpp
    )
    target_include_directories(SPDP_cp_tests PRIVATE tests)
    target_link_libraries(SPDP_cp_tests PRIVATE spdp_core)
    set_target_properties(
        SPDP_cp_tests
        PROPERTIES BUILD_RPATH "${SPDP_BUILD_RPATH}"
    )
    add_test(NAME route_state COMMAND SPDP_cp_tests --case route_state)
    add_test(NAME initial_incumbent_search COMMAND SPDP_cp_tests --case initial_incumbent_search)
    add_test(NAME initial_incumbent_backend COMMAND SPDP_cp_tests --case initial_incumbent_backend)
    add_test(NAME cp_data COMMAND SPDP_cp_tests --case cp_data)
    add_test(NAME cp_log_file COMMAND SPDP_cp_tests --case cp_log_file)
    add_test(NAME cp_symmetry COMMAND SPDP_cp_tests --case cp_symmetry)
    add_test(NAME cp_core COMMAND SPDP_cp_tests --case cp_core)
    add_test(NAME cp_adapter COMMAND SPDP_cp_tests --case cp_adapter)
    add_test(NAME cli COMMAND SPDP_P_step --help)
    set_tests_properties(
        cli
        PROPERTIES PASS_REGULAR_EXPRESSION "initial-incumbent-cp-redundant-terminal-balance"
    )
    add_test(
        NAME cli_invalid_cp_flag
        COMMAND SPDP_P_step --initial-incumbent-cp-redundant-terminal-balance 2
    )
    set_tests_properties(cli_invalid_cp_flag PROPERTIES WILL_FAIL TRUE)

    set(
        SPDP_TEST_DATA_DIR
        "${CMAKE_SOURCE_DIR}/../SPDP_data"
        CACHE PATH
        "Directory containing SPDP benchmark .dat files"
    )
    if(SPDP_ENABLE_LONG_CP_TESTS AND
       EXISTS "${SPDP_TEST_DATA_DIR}/RecDep_day_A6.dat")
        foreach(
            benchmark_case
            benchmark_a6
            benchmark_a6_k2_default
            benchmark_d5_k7
            benchmark_d5_k8
            benchmark_c19_k12
            benchmark_c19_k13
        )
            add_test(
                NAME ${benchmark_case}
                COMMAND SPDP_cp_tests --case ${benchmark_case}
            )
            set_tests_properties(
                ${benchmark_case}
                PROPERTIES
                    LABELS long
                    ENVIRONMENT "SPDP_TEST_DATA_DIR=${SPDP_TEST_DATA_DIR}"
                    TIMEOUT 1300
            )
        endforeach()
    endif()
endif()
