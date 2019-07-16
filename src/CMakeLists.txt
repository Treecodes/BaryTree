# C code for standalone
SET(SRCS_TREEDRIVER_PRIVATE  
    tools.c
    treedriver.c
    partition.c partition.h
    tree-cp.c tree.h
    tree-cp-yuk.c
    tree-pc.c
    tree-pc-coulomb-SS.c
    tree-pc-hermite.c
    tree-pc-hermite-coulombSS.c
    tree-pc-yuk.c
    tree-pc-yuk-SS.c
    tree-batch.c
    xmalloc.c
    sort.c 
    tnode.h
    batch.h
    globvars.h)


SET(SRCS_TREEDRIVER_PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/treedriver.h
    ${CMAKE_CURRENT_SOURCE_DIR}/array.h
    ${CMAKE_CURRENT_SOURCE_DIR}/xmalloc.h
    ${CMAKE_CURRENT_SOURCE_DIR}/tools.h
    ${CMAKE_CURRENT_SOURCE_DIR}/particles.h
    ${CMAKE_CURRENT_SOURCE_DIR}/sort.h)


SET(SRCS_DIRECT
    main-direct.c array.h
    tools.c tools.h 
    xmalloc.c xmalloc.h)


SET(SRCS_DIRECT_DISTRIBUTED
    main-direct-distributed.c array.h
    tools.c tools.h 
    xmalloc.c xmalloc.h)


SET(SRCS_TXT2BIN
    txt2bin.c array.h
    xmalloc.c xmalloc.h
    sort.c)


if(BUILD_TREECODE AND ENABLE_CPU_BUILD)
    add_library(treedriverCPU_ObjLib OBJECT)
    target_sources(treedriverCPU_ObjLib
        PRIVATE ${SRCS_TREEDRIVER_PRIVATE}
        PUBLIC  ${SRCS_TREEDRIVER_PUBLIC})
    set_target_properties(treedriverCPU_ObjLib PROPERTIES
                                               POSITION_INDEPENDENT_CODE ON)
    target_link_libraries(treedriverCPU_ObjLib PUBLIC OpenMP::OpenMP_C
                                               PUBLIC MPI::MPI_C)
    target_compile_features(treedriverCPU_ObjLib PUBLIC c_std_99)


    add_executable(tree-cpu)
    target_sources(tree-cpu PRIVATE main.c)
    target_link_libraries(tree-cpu PRIVATE treedriverCPU_ObjLib)


    add_library(treelib-cpu)
    target_sources(tree-cpu PRIVATE treedriverWrapper.c
                            PUBLIC  treedriverWrapper.h)
    target_link_libraries(treelib-cpu PRIVATE treedriverCPU_ObjLib)


    install(TARGETS tree-cpu treelib-cpu DESTINATION bin)
endif()


if(BUILD_TREECODE AND ENABLE_GPU_BUILD)
#Compute capability all
    add_library(treedriverGPU_ObjLib OBJECT)
    target_sources(treedriverGPU_ObjLib
        PRIVATE ${SRCS_TREEDRIVER_PRIVATE}
        PUBLIC  ${SRCS_TREEDRIVER_PUBLIC})
    set_target_properties(treedriverGPU_ObjLib PROPERTIES
                                               POSITION_INDEPENDENT_CODE ON)
    target_link_libraries(treedriverGPU_ObjLib PUBLIC OpenMP::OpenMP_C
                                               PUBLIC MPI::MPI_C)
    target_compile_features(treedriverGPU_ObjLib PUBLIC c_std_99)
    target_compile_options(treedriverGPU_ObjLib PUBLIC
                                                -acc
                                                -ta=tesla:ccall$<$<CONFIG:Debug>:,time>
                                                $<$<CONFIG:Debug>:-Minfo>)
    target_link_options(treedriverGPU_ObjLib PUBLIC
                                             -acc
                                             -ta=tesla:ccall)


    add_executable(tree-gpu)
    target_sources(tree-gpu PRIVATE main.c)
    target_link_libraries(tree-gpu PRIVATE treedriverGPU_ObjLib)


    add_library(treelib-gpu)
    target_sources(treelib-gpu PRIVATE treedriverWrapper.c
                               PUBLIC  treedriverWrapper.h)
    target_link_libraries(treelib-gpu PRIVATE treedriverGPU_ObjLib)


    install (TARGETS tree-gpu treelib-gpu DESTINATION bin)
endif()


if (BUILD_DIRECT AND ENABLE_CPU_BUILD)
    add_executable(direct-cpu ${SRCS_DIRECT})
    set_target_properties(direct-cpu PROPERTIES
                                     POSITION_INDEPENDENT_CODE ON)
    target_link_libraries(direct-cpu PUBLIC OpenMP::OpenMP_C
                                     PUBLIC MPI::MPI_C)
    target_compile_features(direct-cpu PUBLIC c_std_99)

    add_executable(direct-distributed-cpu ${SRCS_DIRECT_DISTRIBUTED})
    set_target_properties(direct-distributed-cpu PROPERTIES
                                                 POSITION_INDEPENDENT_CODE ON)
    target_link_libraries(direct-distributed-cpu PUBLIC OpenMP::OpenMP_C
                                                 PUBLIC MPI::MPI_C)
    target_compile_features(direct-distributed-cpu PUBLIC c_std_99)

    install (TARGETS direct-cpu direct-distributed-cpu DESTINATION bin)
endif()


if(BUILD_DIRECT AND ENABLE_GPU_BUILD)
#Compute capability all
    add_executable(direct-gpu ${SRCS_DIRECT})
    set_target_properties(direct-gpu PROPERTIES
                                     POSITION_INDEPENDENT_CODE ON)
    target_link_libraries(direct-gpu PUBLIC OpenMP::OpenMP_C
                                     PUBLIC MPI::MPI_C)
    target_compile_features(direct-gpu PUBLIC c_std_99)
    target_compile_options(direct-gpu PUBLIC
                                      -acc
                                      -ta=tesla:ccall$<$<CONFIG:Debug>:,time>
                                      $<$<CONFIG:Debug>:-Minfo>)
    target_link_options(direct-gpu PUBLIC
                                   -acc
                                   -ta=tesla:ccall)

    add_executable(direct-distributed-gpu ${SRCS_DIRECT_DISTRIBUTED})
    set_target_properties(direct-distributed-gpu PROPERTIES
                                                 POSITION_INDEPENDENT_CODE ON)
    target_link_libraries(direct-distributed-gpu PUBLIC OpenMP::OpenMP_C
                                                 PUBLIC MPI::MPI_C)
    target_compile_features(direct-distributed-gpu PUBLIC c_std_99)
    target_compile_options(direct-distributed-gpu PUBLIC
                                                  -acc
                                                  -ta=tesla:ccall$<$<CONFIG:Debug>:,time>
                                                  $<$<CONFIG:Debug>:-Minfo>)
    target_link_options(direct-distributed-gpu PUBLIC
                                               -acc
                                               -ta=tesla:ccall)

    install (TARGETS direct-gpu direct-distributed-gpu DESTINATION bin)
endif()


if (BUILD_TOOLS)
    add_executable(txt2bin ${SRCS_TXT2BIN})
    install (TARGETS txt2bin DESTINATION bin)
endif()
