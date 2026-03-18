cmake_minimum_required(VERSION 3.16)

project(SPDP_P_step_project LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(GUROBI_HOME "/opt/gurobi1301/linux64" CACHE PATH "Path to the Gurobi installation")
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

add_executable(SPDP_P_step
    src/main.cpp
    src/ReadData.cpp
    src/GenMultiGraph.cpp
    src/PstepFormulation.cpp
    ${GUROBI_CPP_SOURCES}
)

target_include_directories(SPDP_P_step PRIVATE include "${GUROBI_INCLUDE_DIR}" "${GUROBI_CPP_DIR}")
target_link_libraries(SPDP_P_step PRIVATE "${GUROBI_CORE_LIB}")

set_target_properties(SPDP_P_step PROPERTIES
    BUILD_RPATH "${GUROBI_LIB_DIR}"
)
