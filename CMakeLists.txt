cmake_minimum_required(VERSION 3.9...3.14)
project(GPU-TREECODE VERSION 0.7 LANGUAGES C)

set(default_build_type "Release")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
        STRING "Choose the type of build." FORCE)

    # Set the possible values of build type for cmake-gui or ccmake
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
                "Debug" "Release" "None")
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)

option(BUILD_TREECODE "build treecode computation executables" ON)
option(BUILD_DIRECT "build direct computation executables" ON)
option(BUILD_TOOLS "build auxiliary utilities" OFF)
option(BUILD_LIBRARIES "build tree libraries for use with external codes" OFF)
option(BUILD_SHARED_LIBS "if building tree libraries, build as shared" ON)

option(ENABLE_GPU_BUILD "build OpenACC GPU version" ON)
option(ENABLE_CPU_BUILD "build CPU version" ON)

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)
set_property(TARGET OpenMP::OpenMP_C APPEND PROPERTY INTERFACE_LINK_LIBRARIES
                                            $<$<C_COMPILER_ID:PGI>:-mp>) 

#Some PGI versions seem to require a link flag for OpenMP
#if(CMAKE_C_COMPILER_ID MATCHES PGI)
#    set(OpenMP_C_LINK_FLAGS -mp)
#endif()

add_subdirectory(src)
