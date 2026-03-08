cmake_minimum_required(VERSION 3.10)

project(SPDP_P_step_project LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(SPDP_P_step
    src/main.cpp
    src/ReadData.cpp
    src/GenMultiGraph.cpp
)

target_include_directories(SPDP_P_step PRIVATE include)
