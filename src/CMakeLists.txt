# C code for standalone
SET(SRCS 
    main.c array.h particles.h tnode.h batch.h globvars.h 
    tools.c tools.h 
    treedriver.c treedriver.h
    partition.c partition.h
    tree-cp.c tree.h
    tree-batch.c
    tree-cp-yuk.c
    tree-pc.c
    tree-pc-hermite.c
    tree-pc-hermite-coulombSS.c
    tree-pc-yuk.c
    tree-pc-yuk-SS.c
    tree-pc-coulomb-SS.c
    xmalloc.c xmalloc.h
    sort.c)

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
    add_executable(tree-cpu ${SRCS})
    set_target_properties(tree-cpu PROPERTIES
                                   POSITION_INDEPENDENT_CODE ON)
    target_link_libraries(tree-cpu PUBLIC OpenMP::OpenMP_C
                                   PUBLIC MPI::MPI_C)
    target_compile_features(tree-cpu PUBLIC c_std_99)

    install (TARGETS tree-cpu DESTINATION bin)
endif()


if(BUILD_TREECODE AND ENABLE_GPU_BUILD)
#Compute capability all
    add_executable(tree-gpu ${SRCS})
    set_target_properties(tree-gpu PROPERTIES
                                   POSITION_INDEPENDENT_CODE ON)
    target_link_libraries(tree-gpu PUBLIC OpenMP::OpenMP_C
                                   PUBLIC MPI::MPI_C)
    target_compile_features(tree-gpu PUBLIC c_std_99)
    target_compile_options(tree-gpu PUBLIC
                                    -acc
                                    -ta=tesla:ccall$<$<CONFIG:Debug>:,time>
                                    $<$<CONFIG:Debug>:-Minfo>)
    target_link_options(tree-gpu PUBLIC
                                 -acc
                                 -ta=tesla:ccall)

    install (TARGETS tree-gpu DESTINATION bin)
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
