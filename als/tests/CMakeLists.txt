# Set the minimum version of CMake that can be used
cmake_minimum_required(VERSION 3.16)

# Set the project name
project(als_tests)

# GoogleTest requires at least C++14
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Download and unpack googletest at configure time
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# Enable testing
enable_testing()

# Add fmt library
find_package(fmt)

# Add a testing executable and link the GoogleTest library
add_executable(
  ${PROJECT_NAME}.out
  als_tests.cc
)
target_link_libraries(
  ${PROJECT_NAME}.out
  PRIVATE gtest_main
  PRIVATE libabc 
  PRIVATE fmt::fmt
  PRIVATE als
)

# Discover tests
include(GoogleTest)
gtest_discover_tests(${PROJECT_NAME}.out)