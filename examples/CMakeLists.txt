if(BUILD_EXAMPLES)

   add_executable(random_cube_cpu random_cube.c)
   target_link_libraries(random_cube_cpu PRIVATE BaryTree_cpu Zoltan_Interface)
   install(TARGETS random_cube_cpu DESTINATION bin)
   
   add_executable(random_cube_large_cpu random_cube_large.c)
   target_link_libraries(random_cube_large_cpu PRIVATE BaryTree_cpu Zoltan_Interface)
   install(TARGETS random_cube_large_cpu DESTINATION bin)

   add_executable(testBaryTreeInterface_cpu testBaryTreeInterface.c)
   target_link_libraries(testBaryTreeInterface_cpu PRIVATE BaryTree_cpu)
   install(TARGETS testBaryTreeInterface_cpu DESTINATION bin)
   
   if(ENABLE_GPU_BUILD)
      add_executable(random_cube_gpu random_cube.c)
      target_link_libraries(random_cube_gpu PRIVATE BaryTree_gpu Zoltan_Interface)
      install(TARGETS random_cube_gpu DESTINATION bin)
   
      add_executable(random_cube_large_gpu random_cube_large.c)
      target_link_libraries(random_cube_large_gpu PRIVATE BaryTree_gpu Zoltan_Interface)
      install(TARGETS random_cube_large_gpu DESTINATION bin)

      add_executable(testBaryTreeInterface_gpu testBaryTreeInterface.c)
      target_link_libraries(testBaryTreeInterface_gpu PRIVATE BaryTree_gpu)
      install(TARGETS testBaryTreeInterface_gpu DESTINATION bin)
   endif()

endif()
