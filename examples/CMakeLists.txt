if(BUILD_EXAMPLES)

   set (AUX_SRCS
        zoltan_fns.h
        zoltan_fns.c
        support_fns.h
        support_fns.c)

   add_executable(random_cube_cpu random_cube.c ${AUX_SRCS})
   target_link_libraries(random_cube_cpu PRIVATE BaryTree_cpu Zoltan_Interface)
   install(TARGETS random_cube_cpu DESTINATION bin)
   
   add_executable(random_cube_reproducible_cpu random_cube_reproducible.c ${AUX_SRCS})
   target_link_libraries(random_cube_reproducible_cpu PRIVATE BaryTree_cpu Zoltan_Interface)
   install(TARGETS random_cube_reproducible_cpu DESTINATION bin)

   add_executable(run_readin_cpu run_readin.c ${AUX_SRCS})
   target_link_libraries(run_readin_cpu PRIVATE BaryTree_cpu Zoltan_Interface)
   install(TARGETS run_readin_cpu DESTINATION bin)

   add_executable(test_BaryTreeInterface_cpu test_BaryTreeInterface.c)
   target_link_libraries(test_BaryTreeInterface_cpu PRIVATE BaryTree_cpu)
   install(TARGETS test_BaryTreeInterface_cpu DESTINATION bin)
   
   if(ENABLE_GPU_BUILD)
      add_executable(random_cube_gpu random_cube.c ${AUX_SRCS})
      target_link_libraries(random_cube_gpu PRIVATE BaryTree_gpu Zoltan_Interface)
      install(TARGETS random_cube_gpu DESTINATION bin)
   
      add_executable(random_cube_reproducible_gpu random_cube_reproducible.c ${AUX_SRCS})
      target_link_libraries(random_cube_reproducible_gpu PRIVATE BaryTree_gpu Zoltan_Interface)
      install(TARGETS random_cube_reproducible_gpu DESTINATION bin)

      add_executable(run_readin_gpu run_readin.c ${AUX_SRCS})
      target_link_libraries(run_readin_gpu PRIVATE BaryTree_gpu Zoltan_Interface)
      install(TARGETS run_readin_gpu DESTINATION bin)

      add_executable(test_BaryTreeInterface_gpu test_BaryTreeInterface.c)
      target_link_libraries(test_BaryTreeInterface_gpu PRIVATE BaryTree_gpu)
      install(TARGETS test_BaryTreeInterface_gpu DESTINATION bin)
   endif()

endif()
