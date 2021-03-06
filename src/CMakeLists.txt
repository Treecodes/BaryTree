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
    utilities/advanced_timings.c
    utilities/advanced_timings.h
# Calling functions for tree and pp calculations
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
# Comm Types struct and related functions
    comm_types/struct_comm_types.h
    comm_types/comm_types.h
    comm_types/comm_types.c
# Comm Windows struct and related functions
    comm_windows/struct_comm_windows.h
    comm_windows/comm_windows.h
    comm_windows/comm_windows.c
# Comm for cluster-particle
    comm_cp/comm_cp.h
    comm_cp/comm_cp.c
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
    interaction_compute/interaction_compute_correction.c
    interaction_compute/interaction_compute_pc.c
    interaction_compute/interaction_compute_cp.c
    interaction_compute/interaction_compute_cc.c
    interaction_compute/interaction_compute_downpass.c)


SET(SRCS_K_COULOMB
# Singularity skipping
    kernels/coulomb/coulomb.h
    kernels/coulomb/coulomb_pp.h
    kernels/coulomb/coulomb_pp.c
    kernels/coulomb/coulomb_pc.h
    kernels/coulomb/coulomb_pc.c
    kernels/coulomb/coulomb_cp.h
    kernels/coulomb/coulomb_cp.c
    kernels/coulomb/coulomb_cc_hermite.h
    kernels/coulomb/coulomb_cc_hermite.c
# Singularity subtraction
    kernels/coulomb/coulomb_ss_pp.h
    kernels/coulomb/coulomb_ss_pp.c
    kernels/coulomb/coulomb_ss_correction.h
    kernels/coulomb/coulomb_ss_correction.c
    kernels/coulomb/coulomb_ss_pc.h
    kernels/coulomb/coulomb_ss_pc.c
    kernels/coulomb/coulomb_ss_cp.h
    kernels/coulomb/coulomb_ss_cp.c
    kernels/coulomb/coulomb_ss_cc.h
    kernels/coulomb/coulomb_ss_cc.c)


SET(SRCS_K_YUKAWA
# Singularity skipping
    kernels/yukawa/yukawa.h
    kernels/yukawa/yukawa_pp.h
    kernels/yukawa/yukawa_pp.c
    kernels/yukawa/yukawa_pc.h
    kernels/yukawa/yukawa_pc.c
    kernels/yukawa/yukawa_cp.h
    kernels/yukawa/yukawa_cp.c
# Singularity subtraction
    kernels/yukawa/yukawa_ss_pp.h
    kernels/yukawa/yukawa_ss_pp.c
    kernels/yukawa/yukawa_ss_correction.h
    kernels/yukawa/yukawa_ss_correction.c
    kernels/yukawa/yukawa_ss_pc.h
    kernels/yukawa/yukawa_ss_pc.c
    kernels/yukawa/yukawa_ss_cp.h
    kernels/yukawa/yukawa_ss_cp.c
    kernels/yukawa/yukawa_ss_cc.h
    kernels/yukawa/yukawa_ss_cc.c)


SET(SRCS_K_REGULARIZED_COULOMB
# Singularity skipping
    kernels/regularized-coulomb/regularized-coulomb.h
    kernels/regularized-coulomb/regularized-coulomb_pp.h
    kernels/regularized-coulomb/regularized-coulomb_pp.c
    kernels/regularized-coulomb/regularized-coulomb_pc.h
    kernels/regularized-coulomb/regularized-coulomb_pc.c
    kernels/regularized-coulomb/regularized-coulomb_cp.h
    kernels/regularized-coulomb/regularized-coulomb_cp.c
# Singularity subtraction
    kernels/regularized-coulomb/regularized-coulomb_ss_pp.h
    kernels/regularized-coulomb/regularized-coulomb_ss_pp.c
    kernels/regularized-coulomb/regularized-coulomb_ss_correction.h
    kernels/regularized-coulomb/regularized-coulomb_ss_correction.c
    kernels/regularized-coulomb/regularized-coulomb_ss_pc.h
    kernels/regularized-coulomb/regularized-coulomb_ss_pc.c)
#   kernels/regularized-coulomb/regularized-coulomb_ss_cp.h
#   kernels/regularized-coulomb/regularized-coulomb_ss_cp.c)


SET(SRCS_K_REGULARIZED_YUKAWA
#Singularity skipping
    kernels/regularized-yukawa/regularized-yukawa.h
    kernels/regularized-yukawa/regularized-yukawa_pp.h
    kernels/regularized-yukawa/regularized-yukawa_pp.c
    kernels/regularized-yukawa/regularized-yukawa_pc.h
    kernels/regularized-yukawa/regularized-yukawa_pc.c
    kernels/regularized-yukawa/regularized-yukawa_cp.h
    kernels/regularized-yukawa/regularized-yukawa_cp.c
#Singularity subtraction
    kernels/regularized-yukawa/regularized-yukawa_ss_pp.h
    kernels/regularized-yukawa/regularized-yukawa_ss_pp.c
    kernels/regularized-yukawa/regularized-yukawa_ss_correction.h
    kernels/regularized-yukawa/regularized-yukawa_ss_correction.c
    kernels/regularized-yukawa/regularized-yukawa_ss_pc.h
    kernels/regularized-yukawa/regularized-yukawa_ss_pc.c)
#   kernels/regularized-yukawa/regularized-yukawa_ss_cp.h
#   kernels/regularized-yukawa/regularized-yukawa_ss_cp.c)


SET(SRCS_K_ATAN
    kernels/atan/atan.h
    kernels/atan/atan_pp.h
    kernels/atan/atan_pp.c
    kernels/atan/atan_pc.h
    kernels/atan/atan_pc.c)
#   kernels/atan/atan_cp.h
#   kernels/atan/atan_cp.c)


SET(SRCS_K_SIN_OVER_R
    kernels/sin-over-r/sin-over-r.h
    kernels/sin-over-r/sin-over-r_pp.h
    kernels/sin-over-r/sin-over-r_pp.c
    kernels/sin-over-r/sin-over-r_pc.h
    kernels/sin-over-r/sin-over-r_pc.c
    kernels/sin-over-r/sin-over-r_cp.h
    kernels/sin-over-r/sin-over-r_cp.c)

SET(SRCS_K_MQ
    kernels/mq/mq.h
    kernels/mq/mq_pp.h
    kernels/mq/mq_pp.c
    kernels/mq/mq_pc.h
    kernels/mq/mq_pc.c)

SET(SRCS_K_RBS_U
    kernels/rbs-u/rbs-u.h
    kernels/rbs-u/rbs-u_pp.h
    kernels/rbs-u/rbs-u_pp.c
    kernels/rbs-u/rbs-u_pc.h
    kernels/rbs-u/rbs-u_pc.c
    kernels/rbs-u/rbs-u_cp.h
    kernels/rbs-u/rbs-u_cp.c)

SET(SRCS_K_RBS_V
    kernels/rbs-v/rbs-v.h
    kernels/rbs-v/rbs-v_pp.h
    kernels/rbs-v/rbs-v_pp.c
    kernels/rbs-v/rbs-v_pc.h
    kernels/rbs-v/rbs-v_pc.c
    kernels/rbs-v/rbs-v_cp.h
    kernels/rbs-v/rbs-v_cp.c)

SET(SRCS_K_USER
    kernels/user_kernel/user_kernel.h
    kernels/user_kernel/user_kernel_pp.h
    kernels/user_kernel/user_kernel_pp.c
    kernels/user_kernel/user_kernel_pc.h
    kernels/user_kernel/user_kernel_pc.c
    kernels/user_kernel/user_kernel_cp.h
    kernels/user_kernel/user_kernel_cp.c)



SET(SRCS_KERNELS ${SRCS_K_COULOMB}
                 ${SRCS_K_YUKAWA}
                 ${SRCS_K_REGULARIZED_COULOMB}
                 ${SRCS_K_REGULARIZED_YUKAWA}
                 ${SRCS_K_ATAN}
                 ${SRCS_K_SIN_OVER_R}
                 ${SRCS_K_MQ}
                 ${SRCS_K_RBS_U}
                 ${SRCS_K_RBS_V}
                 ${SRCS_K_USER})


set(TRGT BaryTree_cpu)
add_library(${TRGT} ${SRCS_TREEDRIVER} ${SRCS_KERNELS})
target_sources(${TRGT}          PRIVATE 
                                    interface/BaryTreeInterface.c
                                PUBLIC  
                                   ${CMAKE_CURRENT_SOURCE_DIR}/interface/BaryTreeInterface.h)
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

install(FILES   ../interfaces/fortran/BaryTreeInterface.fh DESTINATION include)


if(ENABLE_GPU_BUILD)
    set(TRGT BaryTree_gpu)
    add_library(${TRGT} ${SRCS_TREEDRIVER} ${SRCS_KERNELS})
    target_compile_definitions(${TRGT} PRIVATE OPENACC_ENABLED)
    target_sources(${TRGT}        PRIVATE
                                      interface/BaryTreeInterface.c
                                  PUBLIC 
                                      ${CMAKE_CURRENT_SOURCE_DIR}/interface/BaryTreeInterface.h)
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

    install(FILES   ../interfaces/fortran/BaryTreeInterface.fh DESTINATION include)
endif()
