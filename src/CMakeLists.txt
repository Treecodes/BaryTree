SET(SRCS_TREEDRIVER
# Auxiliary utilities
    utilities/timers.h
    utilities/timers.c
    utilities/tools.h
    utilities/tools.c
    utilities/xmalloc.h
    utilities/xmalloc.c
    utilities/array.h
    utilities/enums.h
# Calling functions for tree and direct calculations
    drivers/directdriver.h
    drivers/directdriver.c
    drivers/treedriver.h
    drivers/treedriver.c
# Tree node struct and related functions
    tree/struct_tree.h
    tree/tree.h
    tree/tree.c
    tree/batches.h
    tree/batches.c
    tree/struct_tree_linked_list_node.h
    tree/tree_linked_list.h
    tree/tree_linked_list.c
# Partitioning routines for tree building
    tree/partition.h
    tree/partition.c
# Cluster struct and related functions
    clusters/struct_clusters.h
    clusters/clusters.h
    clusters/clusters.c
# Particle struct and related functions
    particles/struct_particles.h
    particles/particles.h
    particles/particles.c
# Kernel struct and related functions
    run_params/struct_run_params.h
    run_params/run_params.h
    run_params/run_params.c
# Interaction list and related functions
    interaction_lists/struct_interaction_lists.h
    interaction_lists/interaction_lists.h
    interaction_lists/interaction_lists.c
# Interaction compute related functions
    interaction_compute/interaction_compute.h
    interaction_compute/interaction_compute_direct.c
    interaction_compute/interaction_compute_pc.c
    interaction_compute/interaction_compute_cp.c
    interaction_compute/interaction_compute_cc.c
    interaction_compute/interaction_compute_downpass.c)


SET(SRCS_K_COULOMB
# Singularity skipping
    kernels/coulomb/coulomb.h
    kernels/coulomb/coulomb_direct.h
    kernels/coulomb/coulomb_direct.c
    kernels/coulomb/coulomb_pc.h
    kernels/coulomb/coulomb_pc.c
    kernels/coulomb/coulomb_cp.h
    kernels/coulomb/coulomb_cp.c)


SET(SRCS_K_YUKAWA
# Singularity skipping
    kernels/yukawa/yukawa.h
    kernels/yukawa/yukawa_direct.h
    kernels/yukawa/yukawa_direct.c
    kernels/yukawa/yukawa_pc.h
    kernels/yukawa/yukawa_pc.c
    kernels/yukawa/yukawa_cp.h
    kernels/yukawa/yukawa_cp.c)


SET(SRCS_K_TCF
    kernels/tcf/tcf.h
    kernels/tcf/tcf_direct.h
    kernels/tcf/tcf_direct.c
    kernels/tcf/tcf_pc.h
    kernels/tcf/tcf_pc.c
    kernels/tcf/tcf_cp.h
    kernels/tcf/tcf_cp.c)


SET(SRCS_K_DCF
    kernels/dcf/dcf.h
    kernels/dcf/dcf_direct.h
    kernels/dcf/dcf_direct.c
    kernels/dcf/dcf_pc.h
    kernels/dcf/dcf_pc.c
    kernels/dcf/dcf_cp.h
    kernels/dcf/dcf_cp.c)


SET(SRCS_KERNELS ${SRCS_K_COULOMB}
                 ${SRCS_K_YUKAWA}
                 ${SRCS_K_TCF}
                 ${SRCS_K_DCF})


set(TRGT BaryTree_cpu)
add_library(${TRGT} ${SRCS_TREEDRIVER} ${SRCS_KERNELS})
set_target_properties(${TRGT}   PROPERTIES
                                    POSITION_INDEPENDENT_CODE ON)
target_compile_features(${TRGT} PUBLIC
                                    c_std_99)
target_link_libraries(${TRGT}   PUBLIC
                                    OpenMP::OpenMP_C
                                    $<$<C_COMPILER_ID:GNU>:m>)

install(TARGETS ${TRGT}                         LIBRARY DESTINATION lib
                                                ARCHIVE DESTINATION lib)


if(ENABLE_GPU_BUILD)
    set(TRGT BaryTree_gpu)
    add_library(${TRGT} ${SRCS_TREEDRIVER} ${SRCS_KERNELS})
    target_compile_definitions(${TRGT} PRIVATE OPENACC_ENABLED)
    set_target_properties(${TRGT} PROPERTIES
                                      POSITION_INDEPENDENT_CODE ON)
    target_compile_features(${TRGT} PUBLIC
                                      c_std_99)
    target_link_libraries(${TRGT} PUBLIC
                                      OpenMP::OpenMP_C
                                      OpenACC_C
                                      $<$<C_COMPILER_ID:GNU>:m>)

    install(TARGETS ${TRGT}                       LIBRARY DESTINATION lib
                                                  ARCHIVE DESTINATION lib)
endif()
