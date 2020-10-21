add_executable(tests_cpu serial_tests.c)
target_link_libraries(tests_cpu PRIVATE BaryTree_cpu)

add_test(NAME direct_sum_on_10_particles COMMAND tests_cpu 0)
add_test(NAME treecode_on_100_particles COMMAND tests_cpu 1)
add_test(NAME treecode_on_1_target_10000_sources COMMAND tests_cpu 2)
add_test(NAME treecode_parameters_on_1_target_10000_sources COMMAND tests_cpu 3)
add_test(NAME test_treecode_wrapper COMMAND tests_cpu 4)
add_test(NAME test_BLDTT COMMAND tests_cpu 5)