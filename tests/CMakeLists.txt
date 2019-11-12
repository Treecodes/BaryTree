if(BUILD_TESTS AND ENABLE_CPU_BUILD)
   add_executable(tests_cpu tests.c)
   target_link_libraries(tests_cpu PRIVATE treelib-cpu zoltan)
   install(TARGETS tests_cpu DESTINATION bin)
endif()

if(BUILD_TESTS AND ENABLE_GPU_BUILD)
   add_executable(tests_gpu tests.c)
   target_link_libraries(tests_gpu PRIVATE treelib-gpu zoltan)
   install(TARGETS tests_gpu DESTINATION bin)
endif()

if(BUILD_TESTS AND ENABLE_CPU_BUILD)
    set(UNIT_TEST tests_cpu)
    add_custom_command(
         TARGET ${UNIT_TEST}
         COMMENT "Running tests"
         POST_BUILD 
         COMMAND ${UNIT_TEST}
    )
endif()