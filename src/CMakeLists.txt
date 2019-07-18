SET(SRCS_TREEDRIVER
    tools.c
    treedriver.c
    partition.c
    partition.h
    tree-cp.c tree.h
    tree-cp-yuk.c
    tree-pc.c
    tree-pc-hermite.c
    tree-pc-coulomb-SS.c
    tree-pc-hermite-coulomb-SS.c
    tree-pc-yuk.c
    tree-pc-hermite-yuk.c
    tree-pc-yuk-SS.c
    tree-batch.c
    xmalloc.c
    sort.c 
    tnode.h
    batch.h
    globvars.h
    treedriver.h
    array.h
    xmalloc.h
    tools.h
    particles.h
    sort.h)


SET(SRCS_DIRECT
    main-direct.c array.h
    tools.c tools.h 
    xmalloc.c xmalloc.h)


SET(SRCS_DIRECT_DISTRIBUTED
    main-direct-distributed.c array.h
    tools.c tools.h 
    xmalloc.c xmalloc.h)


SET(SRCS_TXT2BIN
    txt2bin.c array.h
    xmalloc.c xmalloc.h
    sort.c)


if(BUILD_TREECODE AND ENABLE_CPU_BUILD)
    set(TRGT tree-cpu)
    add_executable(${TRGT} ${SRCS_TREEDRIVER})
    target_sources(${TRGT}          PRIVATE 
                                        main.c)
    target_compile_features(${TRGT} PRIVATE
                                        c_std_99)
    target_link_libraries(${TRGT}   PRIVATE
                                        OpenMP::OpenMP_C
                                        MPI::MPI_C
                                        $<$<C_COMPILER_ID:GNU>:m>)
    install(TARGETS ${TRGT} DESTINATION bin)
endif()


if(BUILD_TREECODE AND ENABLE_GPU_BUILD)
    set(TRGT tree-gpu)
    add_executable(${TRGT} ${SRCS_TREEDRIVER})
    target_sources(${TRGT}        PRIVATE
                                      main.c)
    target_link_libraries(${TRGT} PRIVATE
                                      OpenMP::OpenMP_C
                                      MPI::MPI_C
                                      OpenACC_C
                                      $<$<C_COMPILER_ID:GNU>:m>)
    install(TARGETS ${TRGT} DESTINATION bin)
endif()


if (BUILD_DIRECT AND ENABLE_CPU_BUILD)
    set(TRGT direct-cpu)
    add_executable(${TRGT} ${SRCS_DIRECT})
    target_compile_features(${TRGT} PRIVATE
                                        c_std_99)
    target_link_libraries(${TRGT}   PRIVATE
                                        OpenMP::OpenMP_C
                                        MPI::MPI_C
                                        $<$<C_COMPILER_ID:GNU>:m>)
    install(TARGETS ${TRGT} DESTINATION bin)


    set(TRGT direct-distributed-cpu)
    add_executable(${TRGT} ${SRCS_DIRECT_DISTRIBUTED})
    target_compile_features(${TRGT} PRIVATE
                                        c_std_99)
    target_link_libraries(${TRGT}   PRIVATE
                                        OpenMP::OpenMP_C
                                        MPI::MPI_C
                                        $<$<C_COMPILER_ID:GNU>:m>)
    install(TARGETS ${TRGT} DESTINATION bin)
endif()


if(BUILD_DIRECT AND ENABLE_GPU_BUILD)
    set(TRGT direct-gpu)
    add_executable(${TRGT} ${SRCS_DIRECT})
    target_link_libraries(${TRGT} PRIVATE
                                      OpenMP::OpenMP_C
                                      MPI::MPI_C
                                      OpenACC_C
                                      $<$<C_COMPILER_ID:GNU>:m>)
    install(TARGETS ${TRGT} DESTINATION bin)


    set(TRGT direct-distributed-gpu)
    add_executable(${TRGT} ${SRCS_DIRECT_DISTRIBUTED})
    target_link_libraries(${TRGT} PRIVATE
                                      OpenMP::OpenMP_C
                                      MPI::MPI_C
                                      OpenACC_C
                                      $<$<C_COMPILER_ID:GNU>:m>)
    install(TARGETS ${TRGT} DESTINATION bin)
endif()


if (BUILD_TOOLS)
    add_executable(txt2bin ${SRCS_TXT2BIN})
    install (TARGETS txt2bin DESTINATION bin)
endif()


if(BUILD_TREE_LIBRARIES AND ENABLE_CPU_BUILD)
    set(TRGT treelib-cpu)
    add_library(${TRGT} ${SRCS_TREEDRIVER})
    target_sources(${TRGT}          PRIVATE 
                                        treedriverWrapper.c
                                     PUBLIC  
                                        ${CMAKE_CURRENT_SOURCE_DIR}/treedriverWrapper.h)
    set_target_properties(${TRGT}   PROPERTIES
                                        POSITION_INDEPENDENT_CODE ON)
    target_compile_features(${TRGT} PUBLIC
                                        c_std_99)
    target_link_libraries(${TRGT}   PUBLIC
                                        OpenMP::OpenMP_C
                                        MPI::MPI_C
                                        $<$<C_COMPILER_ID:GNU>:m>)

    install(TARGETS ${TRGT}               LIBRARY DESTINATION lib)
    install(FILES   treedriverWrapper.h           DESTINATION include)
endif()


if(BUILD_TREE_LIBRARIES AND ENABLE_GPU_BUILD)
    set(TRGT treelib-gpu)
    add_library(${TRGT} ${SRCS_TREEDRIVER})
    target_sources(${TRGT}        PRIVATE
                                      treedriverWrapper.c
                                  PUBLIC 
                                      ${CMAKE_CURRENT_SOURCE_DIR}/treedriverWrapper.h)
    set_target_properties(${TRGT} PROPERTIES
                                      POSITION_INDEPENDENT_CODE ON)
    target_link_libraries(${TRGT} PUBLIC
                                      OpenMP::OpenMP_C
                                      MPI::MPI_C
                                      OpenACC_C
                                      $<$<C_COMPILER_ID:GNU>:m>)

    install(TARGETS ${TRGT}             LIBRARY DESTINATION lib)
    install(FILES   treedriverWrapper.h         DESTINATION include)
endif()
