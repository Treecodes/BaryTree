SET(SRCS_TREEDRIVER
# Auxiliary utilities
    utilities/tools.h
    utilities/tools.c
    utilities/xmalloc.h
    utilities/xmalloc.c
    utilities/array.h
    utilities/enums.h
# Global variables, which will be deprecated
    globvars.h
# Calling functions for tree and direct calculations
    drivers/directdriver.h
    drivers/directdriver.c
    drivers/treedriver.h
    drivers/treedriver.c
# Tree node struct and related functions
    tree/struct_nodes.h
    tree/tree.h
    tree/tree.c
    tree/batches.h
    tree/batches.c
# Partitioning routines for tree building
    tree/partition.h
    tree/partition.c
# Cluster struct and related functions
    clusters/struct_clusters.h
    clusters/clusters.h
    clusters/clusters.c
    clusters/clusters_cp.c
# Particle struct and related functions
    particles/struct_particles.h
    particles/particles.h
    particles/particles.c
# Kernel struct and related functions
    run_params/struct_run_params.h
    run_params/run_params.h
    run_params/run_params.c
# Interaction list related functions
    interaction_lists/interaction_lists.h
    interaction_lists/interaction_lists.c
    interaction_lists/interaction_lists_cp_remote.c
    interaction_lists/interaction_lists_cc.c
# Interaction compute related functions
    interaction_compute/interaction_compute.h
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
    kernels/coulomb/coulomb_cp.c
# Singularity subtraction
    kernels/coulomb/coulomb_ss_direct.h
    kernels/coulomb/coulomb_ss_direct.c
    kernels/coulomb/coulomb_ss_correction.h
    kernels/coulomb/coulomb_ss_correction.c
    kernels/coulomb/coulomb_ss_pc.h
    kernels/coulomb/coulomb_ss_pc.c)
#   kernels/coulomb/coulomb_ss_cp.h
#   kernels/coulomb/coulomb_ss_cp.c)


SET(SRCS_K_YUKAWA
# Singularity skipping
    kernels/yukawa/yukawa.h
    kernels/yukawa/yukawa_direct.h
    kernels/yukawa/yukawa_direct.c
    kernels/yukawa/yukawa_pc.h
    kernels/yukawa/yukawa_pc.c
    kernels/yukawa/yukawa_cp.h
    kernels/yukawa/yukawa_cp.c
# Singularity subtraction
    kernels/yukawa/yukawa_ss_direct.h
    kernels/yukawa/yukawa_ss_direct.c
    kernels/yukawa/yukawa_ss_correction.h
    kernels/yukawa/yukawa_ss_correction.c
    kernels/yukawa/yukawa_ss_pc.h
    kernels/yukawa/yukawa_ss_pc.c)
#   kernels/yukawa/yukawa_ss_cp.h
#   kernels/yukawa/yukawa_ss_cp.c)


SET(SRCS_K_REGULARIZED_COULOMB
# Singularity skipping
    kernels/regularized-coulomb/regularized-coulomb.h
    kernels/regularized-coulomb/regularized-coulomb_direct.h
    kernels/regularized-coulomb/regularized-coulomb_direct.c
    kernels/regularized-coulomb/regularized-coulomb_pc.h
    kernels/regularized-coulomb/regularized-coulomb_pc.c
#   kernels/regularized-coulomb/regularized-coulomb_cp.h
#   kernels/regularized-coulomb/regularized-coulomb_cp.c
# Singularity subtraction
    kernels/regularized-coulomb/regularized-coulomb_ss_direct.h
    kernels/regularized-coulomb/regularized-coulomb_ss_direct.c
    kernels/regularized-coulomb/regularized-coulomb_ss_correction.h
    kernels/regularized-coulomb/regularized-coulomb_ss_correction.c
    kernels/regularized-coulomb/regularized-coulomb_ss_pc.h
    kernels/regularized-coulomb/regularized-coulomb_ss_pc.c)
#   kernels/regularized-coulomb/regularized-coulomb_ss_cp.h
#   kernels/regularized-coulomb/regularized-coulomb_ss_cp.c)


SET(SRCS_K_REGULARIZED_YUKAWA
#Singularity skipping
    kernels/regularized-yukawa/regularized-yukawa.h
    kernels/regularized-yukawa/regularized-yukawa_direct.h
    kernels/regularized-yukawa/regularized-yukawa_direct.c
    kernels/regularized-yukawa/regularized-yukawa_pc.h
    kernels/regularized-yukawa/regularized-yukawa_pc.c
#   kernels/regularized-yukawa/regularized-yukawa_cp.h
#   kernels/regularized-yukawa/regularized-yukawa_cp.c
#Singularity subtraction
    kernels/regularized-yukawa/regularized-yukawa_ss_direct.h
    kernels/regularized-yukawa/regularized-yukawa_ss_direct.c
    kernels/regularized-yukawa/regularized-yukawa_ss_correction.h
    kernels/regularized-yukawa/regularized-yukawa_ss_correction.c
    kernels/regularized-yukawa/regularized-yukawa_ss_pc.h
    kernels/regularized-yukawa/regularized-yukawa_ss_pc.c)
#   kernels/regularized-yukawa/regularized-yukawa_ss_cp.h
#   kernels/regularized-yukawa/regularized-yukawa_ss_cp.c)


SET(SRCS_K_ATAN
    kernels/atan/atan.h
    kernels/atan/atan_direct.h
    kernels/atan/atan_direct.c
    kernels/atan/atan_pc.h
    kernels/atan/atan_pc.c)
#   kernels/atan/atan_cp.h
#   kernels/atan/atan_cp.c)


SET(SRCS_KERNELS ${SRCS_K_COULOMB}
                 ${SRCS_K_YUKAWA}
                 ${SRCS_K_REGULARIZED_COULOMB}
                 ${SRCS_K_REGULARIZED_YUKAWA}
                 ${SRCS_K_ATAN})


set(TRGT BaryTree_cpu)
add_library(${TRGT} ${SRCS_TREEDRIVER} ${SRCS_KERNELS})
target_sources(${TRGT}          PRIVATE 
                                    interface/BaryTreeInterface.c
                                 PUBLIC  
                                    interface/BaryTreeInterface.h)
set_target_properties(${TRGT}   PROPERTIES
                                    POSITION_INDEPENDENT_CODE ON)
target_compile_features(${TRGT} PUBLIC
                                    c_std_99)
target_link_libraries(${TRGT}   PUBLIC
                                    OpenMP::OpenMP_C
                                    MPI::MPI_C
                                    $<$<C_COMPILER_ID:GNU>:m>)

install(TARGETS ${TRGT}                         LIBRARY DESTINATION lib
                                                ARCHIVE DESTINATION lib)
install(FILES   interface/BaryTreeInterface.h           DESTINATION include)


if(ENABLE_GPU_BUILD)
    set(TRGT BaryTree_gpu)
    add_library(${TRGT} ${SRCS_TREEDRIVER} ${SRCS_KERNELS})
    target_compile_definitions(${TRGT} PRIVATE OPENACC_ENABLED)
    target_sources(${TRGT}        PRIVATE
                                      interface/BaryTreeInterface.c
                                  PUBLIC 
                                      interface/BaryTreeInterface.h)
    set_target_properties(${TRGT} PROPERTIES
                                      POSITION_INDEPENDENT_CODE ON)
    target_link_libraries(${TRGT} PUBLIC
                                      OpenMP::OpenMP_C
                                      MPI::MPI_C
                                      OpenACC_C
                                      $<$<C_COMPILER_ID:GNU>:m>)

    install(TARGETS ${TRGT}                       LIBRARY DESTINATION lib
                                                  ARCHIVE DESTINATION lib)
    install(FILES   interface/BaryTreeInterface.h         DESTINATION include)
endif()
