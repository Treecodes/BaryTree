cmake_minimum_required(VERSION 3.9...3.14)
project(BaryTree VERSION 1.0 LANGUAGES C)

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

option(BUILD_TREECODE "build treecode computation executables"                 ON) 
option(BUILD_DIRECT "build direct computation executables"                     ON)
option(BUILD_EXAMPLES "build examples with Zoltan"                             ON)

option(BUILD_TREE_LIBRARIES "build tree libraries for use with external codes" ON)
option(BUILD_SHARED_LIBS "if building tree libraries, build as shared"         ON)

option(ENABLE_GPU_BUILD "build OpenACC GPU version"                            ON)
option(ENABLE_CPU_BUILD "build CPU version"                                    ON)
option(ENABLE_GPU_BUILD_OPENMP "build OpenACC GPU version"                     ON)

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

if(BUILD_EXAMPLES)
    find_package(Zoltan REQUIRED)
endif()

#Some PGI versions seem to require a link flag for OpenMP
set_property(TARGET OpenMP::OpenMP_C APPEND
                                     PROPERTY
                                         INTERFACE_LINK_LIBRARIES
                                         $<$<C_COMPILER_ID:PGI>:-mp>) 

#Creating an interface for OpenACC
add_library(OpenACC_C INTERFACE)
set_property(TARGET OpenACC_C        PROPERTY
                                         INTERFACE_COMPILE_OPTIONS
                                         -acc
                                         -ta=tesla:ccall,cc60$<$<CONFIG:Debug>:,time>
                                         -Minfo)
set_property(TARGET OpenACC_C        PROPERTY
                                         INTERFACE_LINK_LIBRARIES
                                         -acc
                                         -ta=tesla:ccall,cc60)

#target_link_directories(OpenACC_C INTERFACE /sw/ascent/cuda/10.1.105/lib64)
#                                        -lnvToolsExt)
#                                        -lnvToolsExt)

#Creating an interface for OpenMP Acc offload--currently these flags are xl specific
add_library(OpenMP_Acc_C INTERFACE)
set_property(TARGET OpenMP_Acc_C     PROPERTY
                                         INTERFACE_COMPILE_OPTIONS
                                         -qsmp=omp -qoffload
                                         -qinfo)
set_property(TARGET OpenMP_Acc_C     PROPERTY
                                         INTERFACE_LINK_LIBRARIES
                                         -qsmp=omp -qoffload
                                         -qinfo)

if(ENABLE_GPU_BUILD AND NOT (CMAKE_C_COMPILER_ID MATCHES PGI))
    message("GPU build requires PGI compilers! Turning GPU build off.")
    set(ENABLE_GPU_BUILD OFF)
endif()

if(ENABLE_GPU_BUILD_OPENMP AND NOT (CMAKE_C_COMPILER_ID MATCHES XL))
    message("OpenMP GPU build requires XL compilers! Turning OpenMP GPU build off.")
    set(ENABLE_GPU_BUILD_OPENMP OFF)
endif()

add_subdirectory(src)
add_subdirectory(examples)
