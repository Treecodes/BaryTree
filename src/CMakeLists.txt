file(WRITE null.c "")

SET(SRCS_TREEDRIVER_PRIVATE  
    tools.c
    treedriver.c
    partition.c partition.h
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
    globvars.h)


SET(SRCS_TREEDRIVER_PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/treedriver.h
    ${CMAKE_CURRENT_SOURCE_DIR}/array.h
    ${CMAKE_CURRENT_SOURCE_DIR}/xmalloc.h
    ${CMAKE_CURRENT_SOURCE_DIR}/tools.h
    ${CMAKE_CURRENT_SOURCE_DIR}/particles.h
    ${CMAKE_CURRENT_SOURCE_DIR}/sort.h)


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
    set(TRGT treedriverCPU_ObjLib)
    add_library(${TRGT} STATIC null.c)
    target_sources(${TRGT}
        PRIVATE ${SRCS_TREEDRIVER_PRIVATE}
        PUBLIC  ${SRCS_TREEDRIVER_PUBLIC})
#    set_target_properties(${TRGT} PROPERTIES
#                                  POSITION_INDEPENDENT_CODE ON)
    target_link_libraries(${TRGT} PUBLIC OpenMP::OpenMP_C
                                         MPI::MPI_C
                                         $<$<C_COMPILER_ID:GNU>:m>)
    target_compile_features(${TRGT} PUBLIC c_std_99)


    set(TRGT tree-cpu)
    add_executable(${TRGT} null.c)
    target_sources(${TRGT} PRIVATE main.c)
#    set_target_properties(${TRGT} PROPERTIES
#                                  POSITION_INDEPENDENT_CODE ON)
    target_link_libraries(${TRGT} PRIVATE treedriverCPU_ObjLib)
    install(TARGETS ${TRGT} DESTINATION bin)

    if(BUILD_TREE_LIBRARIES)
    set(TRGT treelib-cpu)
    add_library(${TRGT} null.c)
    target_sources(${TRGT} PRIVATE treedriverWrapper.c
                           PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/treedriverWrapper.h)
    set_target_properties(${TRGT} PROPERTIES
                                  POSITION_INDEPENDENT_CODE ON)
    target_link_libraries(${TRGT} PRIVATE treedriverCPU_ObjLib)
    install(TARGETS ${TRGT}
            LIBRARY DESTINATION lib)
    install(FILES treedriverWrapper.h 
            DESTINATION include)
    endif()
endif()


if(BUILD_TREECODE AND ENABLE_GPU_BUILD)
#Compute capability all
    set(TRGT treedriverGPU_ObjLib)
    add_library(${TRGT} STATIC null.c)
    target_sources(${TRGT}
        PRIVATE ${SRCS_TREEDRIVER_PRIVATE}
        PUBLIC  ${SRCS_TREEDRIVER_PUBLIC})
    set_target_properties(${TRGT} PROPERTIES
#                                 POSITION_INDEPENDENT_CODE ON
                                  LINK_FLAGS "-acc -ta=tesla:ccall,cc60")
    target_link_libraries(${TRGT} PUBLIC OpenMP::OpenMP_C
                                         MPI::MPI_C
                                         $<$<C_COMPILER_ID:GNU>:m>)
    target_compile_options(${TRGT} PUBLIC
                                   -acc
                                   -ta=tesla:ccall,cc60$<$<CONFIG:Debug>:,time>
                                   $<$<CONFIG:Debug>:-Minfo>)
    #target_link_options(${TRGT} PUBLIC -acc -ta=tesla:ccall,cc60)


    set(TRGT tree-gpu)
    add_executable(${TRGT} null.c)
    target_sources(${TRGT} PRIVATE main.c)
    set_target_properties(${TRGT} PROPERTIES
#                                  POSITION_INDEPENDENT_CODE ON
                                  LINK_FLAGS "-acc -ta=tesla:ccall,cc60")
    target_link_libraries(${TRGT} PRIVATE treedriverGPU_ObjLib)
    install(TARGETS ${TRGT} DESTINATION bin)


    if(BUILD_TREE_LIBRARIES)
    set(TRGT treelib-gpu)
    add_library(${TRGT} null.c)
    target_sources(${TRGT} PRIVATE treedriverWrapper.c
                           PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/treedriverWrapper.h)
    set_target_properties(${TRGT} PROPERTIES
                                  POSITION_INDEPENDENT_CODE ON
                                  LINK_FLAGS "-acc -ta=tesla:ccall,cc60")
    target_link_libraries(${TRGT} PRIVATE treedriverGPU_ObjLib)
    install(TARGETS ${TRGT}
            LIBRARY DESTINATION lib)
    install(FILES treedriverWrapper.h 
            DESTINATION include)
    endif()
endif()


if (BUILD_DIRECT AND ENABLE_CPU_BUILD)
    set(TRGT direct-cpu)
    add_executable(${TRGT} ${SRCS_DIRECT})
#    set_target_properties(${TRGT} PROPERTIES
#                                  POSITION_INDEPENDENT_CODE ON)
    target_link_libraries(${TRGT} PRIVATE OpenMP::OpenMP_C
                                          MPI::MPI_C
                                          $<$<C_COMPILER_ID:GNU>:m>)
    target_compile_features(${TRGT} PRIVATE c_std_99)
    install(TARGETS ${TRGT} DESTINATION bin)


    set(TRGT direct-distributed-cpu)
    add_executable(${TRGT} ${SRCS_DIRECT_DISTRIBUTED})
#    set_target_properties(${TRGT} PROPERTIES
#                                  POSITION_INDEPENDENT_CODE ON)
    target_link_libraries(${TRGT} PRIVATE OpenMP::OpenMP_C
                                          MPI::MPI_C
                                          $<$<C_COMPILER_ID:GNU>:m>)
    target_compile_features(${TRGT} PRIVATE c_std_99)
    install(TARGETS ${TRGT} DESTINATION bin)
endif()


if(BUILD_DIRECT AND ENABLE_GPU_BUILD)
#Compute capability all
    set(TRGT direct-gpu)
    add_executable(${TRGT} ${SRCS_DIRECT})
    set_target_properties(${TRGT} PROPERTIES
#                                  POSITION_INDEPENDENT_CODE ON
                                  LINK_FLAGS "-acc -ta=tesla:ccall,cc60")
    target_link_libraries(${TRGT} PRIVATE OpenMP::OpenMP_C
                                          MPI::MPI_C
                                          $<$<C_COMPILER_ID:GNU>:m>)
    target_compile_options(${TRGT} PRIVATE
                                   -acc
                                   -ta=tesla:ccall,cc60$<$<CONFIG:Debug>:,time>
                                   $<$<CONFIG:Debug>:-Minfo>)
    #target_link_options(${TRGT} PRIVATE -acc -ta=tesla:ccall,cc60)
    install(TARGETS ${TRGT} DESTINATION bin)


    set(TRGT direct-distributed-gpu)
    add_executable(${TRGT} ${SRCS_DIRECT_DISTRIBUTED})
    set_target_properties(${TRGT} PROPERTIES
#                                  POSITION_INDEPENDENT_CODE ON
                                  LINK_FLAGS "-acc -ta=tesla:ccall,cc60")
    target_link_libraries(${TRGT} PRIVATE OpenMP::OpenMP_C
                                          MPI::MPI_C
                                          $<$<C_COMPILER_ID:GNU>:m>)
    target_compile_options(${TRGT} PRIVATE
                                   -acc
                                   -ta=tesla:ccall,cc60$<$<CONFIG:Debug>:,time>
                                   $<$<CONFIG:Debug>:-Minfo>)
    #target_link_options(${TRGT} PRIVATE -acc -ta=tesla:ccall,cc60)
    install(TARGETS ${TRGT} DESTINATION bin)
endif()


if (BUILD_TOOLS)
    add_executable(txt2bin ${SRCS_TXT2BIN})
    install (TARGETS txt2bin DESTINATION bin)
endif()
