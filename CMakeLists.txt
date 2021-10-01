cmake_minimum_required(VERSION 3.18)
project(BaryTree VERSION 1.0 LANGUAGES C CXX CUDA)

set(default_build_type "Release")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
        STRING "Choose the type of build." FORCE)

    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
                "Debug" "Release" "None")
endif()
 
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)

option(BUILD_EXAMPLES "build examples"                          ON)
option(ENABLE_GPU_BUILD "build OpenACC GPU version"             ON)
option(ENABLE_CUDA_BUILD "build OpenACC+CUDA GPU version"       ON)
option(BUILD_SHARED_LIBS "build libraries as shared objects"	ON)

find_package(OpenMP REQUIRED)

#Creating an interface for OpenMP
#Some PGI versions seem to require a link flag for OpenMP
set_property(TARGET OpenMP::OpenMP_C APPEND
                                     PROPERTY
                                         INTERFACE_LINK_LIBRARIES
                                         $<$<C_COMPILER_ID:PGI>:-mp>) 

#CUDA build requires GPU build
if(ENABLE_CUDA_BUILD AND NOT ENABLE_GPU_BUILD)
    message(STATUS "CUDA build requires GPU build! Turning CUDA build off.")
    set(ENABLE_CUDA_BUILD OFF CACHE BOOL "build OpenACC+CUDA GPU version" FORCE)
endif()

#GPU build requires PGI for C compiler
#Note that both pgcc and nvc report PGI
if(ENABLE_GPU_BUILD AND NOT (CMAKE_C_COMPILER_ID MATCHES PGI))
    message(STATUS "GPU build requires PGI compilers! Turning GPU build off.")
    set(ENABLE_GPU_BUILD OFF CACHE BOOL "build OpenACC GPU version" FORCE)
endif()

if(ENABLE_CUDA_BUILD)
    enable_language(CUDA)
    set(CMAKE_CUDA_HOST_LINK_LAUNCHER ${CMAKE_CXX_COMPILER})
endif()

if(ENABLE_GPU_BUILD)
    #Creating an interface for OpenACC
    add_library(OpenACC_C INTERFACE)
    set_property(TARGET OpenACC_C PROPERTY
                                      INTERFACE_COMPILE_OPTIONS
                                      $<$<NOT:$<COMPILE_LANGUAGE:CUDA>>:
                                      -acc
                                      -ta=tesla:cc70$<$<CONFIG:Debug>:,time>
                                      -Minfo>)
    set_property(TARGET OpenACC_C PROPERTY
                                      INTERFACE_LINK_LIBRARIES
                                      -acc
                                      -ta=tesla:cc70)
endif()

add_subdirectory(src)
add_subdirectory(examples)
