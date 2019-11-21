if(BUILD_TESTS AND ENABLE_CPU_BUILD)
   add_executable(tests_cpu serial_tests.c)
   target_link_libraries(tests_cpu PRIVATE treelib-cpu Zoltan_Interface)
   install(TARGETS tests_cpu DESTINATION bin)
endif()

if(BUILD_TESTS AND ENABLE_CPU_BUILD)
   add_executable(tests_cpu_mpi mpi_tests.c)
   target_link_libraries(tests_cpu_mpi PRIVATE treelib-cpu Zoltan_Interface)
   install(TARGETS tests_cpu_mpi DESTINATION bin)
endif()

if(BUILD_TESTS AND ENABLE_GPU_BUILD)
   add_executable(tests_gpu serial_tests.c)
   target_link_libraries(tests_gpu PRIVATE treelib-gpu Zoltan_Interface)
   install(TARGETS tests_gpu DESTINATION bin)
endif()

if(BUILD_TESTS AND RUN_TESTS_ON_BUILD AND ENABLE_CPU_BUILD)
   set(UNIT_TEST tests_cpu)
   add_custom_command(
         TARGET ${UNIT_TEST}
         COMMENT "Running tests on CPU"
         POST_BUILD 
         COMMAND ${UNIT_TEST}
    )
endif()

if(BUILD_TESTS AND RUN_TESTS_ON_BUILD AND ENABLE_CPU_BUILD)
   set(UNIT_TEST_MPI tests_cpu_mpi)
   add_custom_command(
         TARGET ${UNIT_TEST_MPI}
         COMMENT "Running tests on CPU"
         POST_BUILD 
         COMMAND mpirun -np 2 ${UNIT_TEST}
    )
endif()
