cmake_minimum_required(VERSION 3.16)

project(SPDP_P_step_project LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(GUROBI_HOME "$ENV{GUROBI_HOME}" CACHE PATH "Path to the Gurobi installation, e.g. /opt/gurobi1301/linux64")

if(NOT GUROBI_HOME)
    message(FATAL_ERROR "GUROBI_HOME is not set. Example: export GUROBI_HOME=/opt/gurobi1301/linux64")
endif()

set(GUROBI_INCLUDE_DIR "${GUROBI_HOME}/include")
set(GUROBI_CPP_DIR "${GUROBI_HOME}/src/cpp")

if(NOT EXISTS "${GUROBI_INCLUDE_DIR}/gurobi_c++.h")
    message(FATAL_ERROR "Could not find gurobi_c++.h under ${GUROBI_INCLUDE_DIR}")
endif()

if(NOT EXISTS "${GUROBI_CPP_DIR}")
    message(FATAL_ERROR "Could not find the Gurobi C++ source directory under ${GUROBI_CPP_DIR}")
endif()

file(GLOB GUROBI_CPP_SOURCES "${GUROBI_CPP_DIR}/*.cpp")

if(NOT GUROBI_CPP_SOURCES)
    message(FATAL_ERROR "Could not find Gurobi C++ wrapper sources under ${GUROBI_CPP_DIR}")
endif()

file(GLOB GUROBI_VERSIONED_LINKS "${GUROBI_HOME}/lib/libgurobi[0-9]*.so")

if(NOT GUROBI_VERSIONED_LINKS)
    message(FATAL_ERROR "Could not find a versioned Gurobi shared library link under ${GUROBI_HOME}/lib")
endif()

list(SORT GUROBI_VERSIONED_LINKS COMPARE NATURAL ORDER DESCENDING)
list(GET GUROBI_VERSIONED_LINKS 0 GUROBI_CORE_LIB)

add_executable(SPDP_P_step
    src/main.cpp
    src/ReadData.cpp
    src/GenMultiGraph.cpp
    ${GUROBI_CPP_SOURCES}
)

target_include_directories(SPDP_P_step PRIVATE include "${GUROBI_INCLUDE_DIR}" "${GUROBI_CPP_DIR}")
target_link_libraries(SPDP_P_step PRIVATE "${GUROBI_CORE_LIB}")

set_target_properties(SPDP_P_step PROPERTIES
    BUILD_RPATH "${GUROBI_HOME}/lib"
)
