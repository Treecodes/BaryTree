SET(SRCS_TREEDRIVER
    tools.c
    tools.h
    xmalloc.c
    xmalloc.h
    array.h
    globvars.h
    directdriver.c
    directdriver.h
    treedriver.c
    treedriver.h
    tree.c
    tree.h
    batches.c
    batches.h
    clusters.c
    clusters.h
    particles.c
    particles.h
    interaction_lists.c
    interaction_lists.h
    interaction_compute.c
    interaction_compute.h
    partition.c
    partition.h
    struct_particles.h
    struct_clusters.h
    struct_nodes.h
    kernels/coulomb.c
    kernels/coulomb.h
    kernels/coulomb_singularity_subtraction.c
    kernels/coulomb_singularity_subtraction.h
    kernels/yukawa.c
    kernels/yukawa.h
    kernels/yukawa_singularity_subtraction.c
    kernels/yukawa_singularity_subtraction.h
    kernels/dcf.c
    kernels/dcf.h
    kernels/tcf.c
    kernels/tcf.h)


SET(SRCS_DIRECT_DISTRIBUTED
    main_direct.c
    array.h
    tools.c
    tools.h 
    xmalloc.c
    xmalloc.h
    kernels/kernels_direct.h
    kernels/kernels_direct.c)


if(BUILD_TREECODE AND ENABLE_CPU_BUILD)
    set(TRGT tree-distributed-cpu)
    add_executable(${TRGT} ${SRCS_TREEDRIVER})
    target_sources(${TRGT}          PRIVATE 
                                        main.c)
    target_compile_features(${TRGT} PRIVATE
                                        c_std_99)
    target_link_libraries(${TRGT}   PRIVATE
                                        OpenMP::OpenMP_C
                                        MPI::MPI_C
                                        $<$<C_COMPILER_ID:GNU>:m>)
    install(TARGETS ${TRGT} DESTINATION bin)
endif()


if(BUILD_TREECODE AND ENABLE_GPU_BUILD)
    set(TRGT tree-distributed-gpu)
    add_executable(${TRGT} ${SRCS_TREEDRIVER})
    target_compile_definitions(${TRGT} PRIVATE OPENACC_ENABLED)
    target_sources(${TRGT}        PRIVATE
                                      main.c)
    target_link_libraries(${TRGT} PRIVATE
                                      OpenMP::OpenMP_C
                                      MPI::MPI_C
                                      OpenACC_C
                                      $<$<C_COMPILER_ID:GNU>:m>)
    install(TARGETS ${TRGT} DESTINATION bin)
endif()


if (BUILD_DIRECT AND ENABLE_CPU_BUILD)
    set(TRGT direct-distributed-cpu)
    add_executable(${TRGT} ${SRCS_DIRECT_DISTRIBUTED})
    target_compile_features(${TRGT} PRIVATE
                                        c_std_99)
    target_link_libraries(${TRGT}   PRIVATE
                                        OpenMP::OpenMP_C
                                        MPI::MPI_C
                                        $<$<C_COMPILER_ID:GNU>:m>)
    install(TARGETS ${TRGT} DESTINATION bin)
endif()


if(BUILD_DIRECT AND ENABLE_GPU_BUILD)
    set(TRGT direct-distributed-gpu)
    add_executable(${TRGT} ${SRCS_DIRECT_DISTRIBUTED})
    target_compile_definitions(${TRGT} PRIVATE OPENACC_ENABLED)
    target_link_libraries(${TRGT} PRIVATE
                                      OpenMP::OpenMP_C
                                      MPI::MPI_C
                                      OpenACC_C
                                      $<$<C_COMPILER_ID:GNU>:m>)
    install(TARGETS ${TRGT} DESTINATION bin)
endif()


if(BUILD_TREE_LIBRARIES AND ENABLE_CPU_BUILD)
    set(TRGT treelib-cpu)
    add_library(${TRGT} ${SRCS_TREEDRIVER})
    target_sources(${TRGT}          PRIVATE 
                                        treedriverWrapper.c
                                     PUBLIC  
                                        ${CMAKE_CURRENT_SOURCE_DIR}/treedriverWrapper.h)
    set_target_properties(${TRGT}   PROPERTIES
                                        POSITION_INDEPENDENT_CODE ON)
    target_compile_features(${TRGT} PUBLIC
                                        c_std_99)
    target_link_libraries(${TRGT}   PUBLIC
                                        OpenMP::OpenMP_C
                                        MPI::MPI_C
                                        $<$<C_COMPILER_ID:GNU>:m>)

    install(TARGETS ${TRGT}               LIBRARY DESTINATION lib
                                          ARCHIVE DESTINATION lib)
    install(FILES   treedriverWrapper.h           DESTINATION include)
endif()


if(BUILD_TREE_LIBRARIES AND ENABLE_GPU_BUILD)
    set(TRGT treelib-gpu)
    add_library(${TRGT} ${SRCS_TREEDRIVER})
    target_compile_definitions(${TRGT} PRIVATE OPENACC_ENABLED)
    target_sources(${TRGT}        PRIVATE
                                      treedriverWrapper.c
                                  PUBLIC 
                                      ${CMAKE_CURRENT_SOURCE_DIR}/treedriverWrapper.h)
    set_target_properties(${TRGT} PROPERTIES
                                      POSITION_INDEPENDENT_CODE ON)
    target_link_libraries(${TRGT} PUBLIC
                                      OpenMP::OpenMP_C
                                      MPI::MPI_C
                                      OpenACC_C
                                      $<$<C_COMPILER_ID:GNU>:m>)

    install(TARGETS ${TRGT}             LIBRARY DESTINATION lib
                                        ARCHIVE DESTINATION lib)
    install(FILES   treedriverWrapper.h         DESTINATION include)
endif()
