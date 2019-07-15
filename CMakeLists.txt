cmake_minimum_required(VERSION 3.1...3.14)
project(GPU-TREECODE VERSION 0.0 LANGUAGES C)

set(default_build_type "Release")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
        STRING "Choose the type of build." FORCE)

    # Set the possible values of build type for cmake-gui or ccmake
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
                "Debug" "Release")
endif()

set(PROJECT_ROOT ${PROJECT_SOURCE_DIR})
set(PROJECT_BUILD ${PROJECT_BINARY_DIR})
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BUILD}/bin)
set(EXECUTABLE_INSTALL_PATH bin)

option(ENABLE_GPU_BUILD "build OpenACC GPU version" ON)
option(ENABLE_CPU_BUILD "build CPU version" ON)
file(MAKE_DIRECTORY build/bin)

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

if(ENABLE_GPU_BUILD)
    find_package(OpenACC REQUIRED)
endif()

add_subdirectory(src)
