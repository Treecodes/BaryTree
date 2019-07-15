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

if(ENABLE_CPU_BUILD)
    add_executable(tree-cpu ${SRCS})
    set_target_properties(tree-cpu PROPERTIES
                                   POSITION_INDEPENDENT_CODE ON)
    target_link_libraries(tree-cpu PUBLIC OpenMP::OpenMP_C
                                   PUBLIC MPI::MPI_C)
    target_compile_features(tree-cpu PUBLIC c_std_99)

    install (TARGETS tree-cpu DESTINATION bin)
endif()

if(ENABLE_GPU_BUILD)
# Compute capability 3.5
    add_executable(tree-gpu-35 ${SRCS})
    set_target_properties(tree-gpu-35 PROPERTIES
                                      POSITION_INDEPENDENT_CODE ON)
    target_link_libraries(tree-gpu-35 PUBLIC OpenMP::OpenMP_C
                                      PUBLIC MPI::MPI_C)
    target_compile_features(tree-gpu-35 PUBLIC c_std_99)
    target_compile_options(tree-gpu-35 PUBLIC $<$<CONFIG:Debug>:-Minfo>)
    target_link_options(tree-gpu-35 PUBLIC -acc -ta=tesla:cc35)

#Compute capability 7.0
    add_executable(tree-gpu-70 ${SRCS})
    set_target_properties(tree-gpu-70 PROPERTIES
                                      POSITION_INDEPENDENT_CODE ON)
    target_link_libraries(tree-gpu-70 PUBLIC OpenMP::OpenMP_C
                                      PUBLIC MPI::MPI_C)
    target_compile_features(tree-gpu-70 PUBLIC c_std_99)
    target_compile_options(tree-gpu-70 PUBLIC $<$<CONFIG:Debug>:-Minfo>)
    target_link_options(tree-gpu-70 PUBLIC -acc -ta=tesla:cc70)

    install (TARGETS tree-gpu-35 tree-gpu-70 DESTINATION bin)
endif()
