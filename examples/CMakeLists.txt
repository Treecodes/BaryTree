if(BUILD_EXAMPLES)

   set (AUX_SRCS
        support_fns.h
        support_fns.c)

   add_executable(run_readin_rism_cpu run_readin_rism.c ${AUX_SRCS})
   target_link_libraries(run_readin_rism_cpu PRIVATE BaryTree_cpu)
   install(TARGETS run_readin_rism_cpu DESTINATION bin)

   if(ENABLE_GPU_BUILD)
      add_executable(run_readin_rism_gpu run_readin_rism.c ${AUX_SRCS})
      target_link_libraries(run_readin_rism_gpu PRIVATE BaryTree_gpu)
      install(TARGETS run_readin_rism_gpu DESTINATION bin)
   endif()

endif()
