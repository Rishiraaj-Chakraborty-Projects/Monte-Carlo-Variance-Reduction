cmake_minimum_required(VERSION 3.15)
project(mc LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_executable(mc src/main.cpp)

if(MSVC)
  target_compile_options(mc PRIVATE /W4 /O2)
else()
  target_compile_options(mc PRIVATE -Wall -Wextra -O2)
endif()


enable_testing()
add_test(NAME selftest COMMAND mc --test)
