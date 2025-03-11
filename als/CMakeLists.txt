# Set the minimum version of CMake that can be used
cmake_minimum_required(VERSION 3.16)

# Set the project name
project(als)

# Find fmt library, if not found, download it
# find_package(fmt QUIET)
# if (NOT fmt_FOUND)
    # message(STATUS "fmt library not found, downloading it")
    include(FetchContent)
    FetchContent_Declare(
        fmt
        GIT_REPOSITORY https://github.com/fmtlib/fmt
        GIT_TAG        e69e5f977d458f2650bb346dadf2ad30c5320281) # 10.2.1
    FetchContent_MakeAvailable(fmt)
# endif()

# Collect source files
aux_source_directory(${PROJECT_SOURCE_DIR}/src sources)

# Add include paths
include_directories(${TOP_DIR}/abc/src)
include_directories(${TOP_DIR}/cryptominisat/src)
include_directories(${PROJECT_SOURCE_DIR}/src)

# Set C++20 standard for the project
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set output path
set(EXECUTABLE_OUTPUT_PATH ${TOP_DIR})
message(STATUS "EXECUTABLE_OUTPUT_PATH = ${EXECUTABLE_OUTPUT_PATH}")

# Create an object library from the source files
add_library(als_objs OBJECT ${sources})
target_compile_features(als_objs PUBLIC cxx_std_20)
target_compile_options(als_objs PRIVATE -Wall -O3 -g)
target_include_directories(als_objs PRIVATE $<TARGET_PROPERTY:fmt::fmt,INTERFACE_INCLUDE_DIRECTORIES>)
target_compile_definitions(als_objs PRIVATE LIN64) # for abc

# Create an executable by linking the object files from the object library
add_executable(${PROJECT_NAME}.out $<TARGET_OBJECTS:als_objs>)
target_link_libraries(${PROJECT_NAME}.out 
    PRIVATE libabc 
    PRIVATE fmt::fmt
    PRIVATE ncurses
    PRIVATE cryptominisat5
    PRIVATE tcmalloc_and_profiler
)

# Create a static library by linking the same object files
add_library(${PROJECT_NAME} STATIC $<TARGET_OBJECTS:als_objs>)
target_link_libraries(${PROJECT_NAME} 
    PRIVATE libabc 
    PRIVATE fmt::fmt
    PRIVATE ncurses
    PRIVATE cryptominisat5
)

# google test
message(STATUS "BUILD_TESTING = ${BUILD_TESTING}")
if (BUILD_TESTING)
    add_subdirectory(tests)
else()
    message(STATUS "Testing is disabled")
endif()