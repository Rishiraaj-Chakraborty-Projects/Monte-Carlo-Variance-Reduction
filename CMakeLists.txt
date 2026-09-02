cmake_minimum_required(VERSION 3.15)
project(MonteCarloVarianceReduction LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(MSVC)
    add_compile_options(/W4 /permissive-)
else()
    add_compile_options(-Wall -Wextra -Wpedantic -Wconversion -Wshadow)
endif()

# Core Library Target
add_library(mc_lib src/simulator.cpp)
target_include_directories(mc_lib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# CLI Demo Target
add_executable(mc_demo apps/main.cpp)
target_link_libraries(mc_demo PRIVATE mc_lib)

# Unit Test Target
enable_testing()
add_executable(mc_tests tests/test_monte_carlo.cpp)
target_link_libraries(mc_tests PRIVATE mc_lib)
add_test(NAME MonteCarloUnitTest COMMAND mc_tests)

# Benchmark Target
add_executable(mc_benchmark benchmarks/bench_monte_carlo.cpp)
target_link_libraries(mc_benchmark PRIVATE mc_lib)
