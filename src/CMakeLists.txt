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
    kernel.c
    kernel.h
    struct_kernel.h
    interaction_lists.c
    interaction_lists.h
    interaction_compute_pc.c
    interaction_compute.h
    partition.c
    partition.h
    struct_particles.h
    struct_clusters.h
    struct_nodes.h
    kernels/coulomb.c
    kernels/coulomb.h
    kernels/regularized-coulomb.c
    kernels/regularized-coulomb.h
    kernels/coulomb_singularity_subtraction.c
    kernels/coulomb_singularity_subtraction.h
    kernels/regularized-coulomb-singularity-subtraction.c
    kernels/regularized-coulomb-singularity-subtraction.h
    kernels/yukawa.c
    kernels/yukawa.h
    kernels/regularized-yukawa.c
    kernels/regularized-yukawa.h
    kernels/yukawa_singularity_subtraction.c
    kernels/yukawa_singularity_subtraction.h
    kernels/regularized-yukawa-singularity-subtraction.c
    kernels/regularized-yukawa-singularity-subtraction.h
    kernels/dcf.c
    kernels/dcf.h
    kernels/tcf.c
    kernels/tcf.h
    kernels/atan.c
    kernels/atan.h
    clusters_cp.c
    interaction_lists_cp.c
    interaction_compute_cp.c
    interaction_lists_cc.c
    interaction_compute_cc.c
    kernels/cp_coulomb.c
    kernels/cp_coulomb.h)


set(TRGT BaryTree_cpu)
add_library(${TRGT} ${SRCS_TREEDRIVER})
target_sources(${TRGT}          PRIVATE 
                                    BaryTreeInterface.c
                                 PUBLIC  
                                    ${CMAKE_CURRENT_SOURCE_DIR}/BaryTreeInterface.h)
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
install(FILES   BaryTreeInterface.h           DESTINATION include)


if(ENABLE_GPU_BUILD)
    set(TRGT BaryTree_gpu)
    add_library(${TRGT} ${SRCS_TREEDRIVER})
    target_compile_definitions(${TRGT} PRIVATE OPENACC_ENABLED)
    target_sources(${TRGT}        PRIVATE
                                      BaryTreeInterface.c
                                  PUBLIC 
                                      ${CMAKE_CURRENT_SOURCE_DIR}/BaryTreeInterface.h)
    set_target_properties(${TRGT} PROPERTIES
                                      POSITION_INDEPENDENT_CODE ON)
    target_link_libraries(${TRGT} PUBLIC
                                      OpenMP::OpenMP_C
                                      MPI::MPI_C
                                      OpenACC_C
                                      $<$<C_COMPILER_ID:GNU>:m>)

    install(TARGETS ${TRGT}             LIBRARY DESTINATION lib
                                        ARCHIVE DESTINATION lib)
    install(FILES   BaryTreeInterface.h         DESTINATION include)
endif()
