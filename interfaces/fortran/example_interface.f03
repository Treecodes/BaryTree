PROGRAM BaryTree_Fortran_Example

    USE MPI
    USE ISO_C_BINDING, ONLY: C_LOC
    IMPLICIT NONE

    INCLUDE "BaryTreeInterface.fh"

    INTEGER :: rank, num_proc, ierr

    INTEGER :: num_targets, num_sources, kernel, num_kernel_params, singularity, &
            approximation, compute_type, degree, max_source_leaf, max_target_leaf, & 
            verbosity

    DOUBLE PRECISION, POINTER, DIMENSION(:) :: target_x, target_y, target_z, target_q, &
            source_x, source_y, source_z, source_q, source_w, potential

    DOUBLE PRECISION, DIMENSION(2) :: kernel_params

    DOUBLE PRECISION :: theta, size_check, beta

    CALL MPI_INIT(ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, num_proc, ierr)

    kernel = RBS_U
    num_kernel_params = 1
    kernel_params(1) = 0.5

    singularity = SKIPPING
    approximation = LAGRANGE
    verbosity = 3

    max_source_leaf = 50
    max_target_leaf = 50
    size_check = 1.0
    beta = -1

    compute_type = CLUSTER_CLUSTER
    theta = 0.7
    degree = 3

    num_targets = 10000
    num_sources = 10000

    ALLOCATE(target_x(num_targets), target_y(num_targets), target_z(num_targets), &
             target_q(num_targets), potential(num_targets), &
             source_x(num_sources), source_y(num_sources), source_z(num_sources), &
             source_q(num_sources), source_w(num_sources))

    CALL RANDOM_NUMBER(target_x)
    CALL RANDOM_NUMBER(target_y)
    CALL RANDOM_NUMBER(target_z)
    CALL RANDOM_NUMBER(target_q)

    CALL RANDOM_NUMBER(source_x)
    CALL RANDOM_NUMBER(source_y)
    CALL RANDOM_NUMBER(source_z)
    CALL RANDOM_NUMBER(source_q)
    CALL RANDOM_NUMBER(source_w)

    ! Calling with kernel as RBS_U
    CALL BaryTreeInterface(num_targets, num_sources, &
             C_LOC(target_x), C_LOC(target_y), C_LOC(target_z), C_LOC(target_q), &
             C_LOC(source_x), C_LOC(source_y), C_LOC(source_z), C_LOC(source_q), &
             C_LOC(source_w), C_LOC(potential), &
             kernel, num_kernel_params, C_LOC(kernel_params), &
             singularity, approximation, compute_type, theta, degree, &
             max_source_leaf, max_target_leaf, size_check, beta, &
             verbosity)

    PRINT *, "RBS u total potential is: ", SUM(potential)

    ! Calling with kernel as RBS_V
    CALL BaryTreeInterface(num_targets, num_sources, &
             C_LOC(target_x), C_LOC(target_y), C_LOC(target_z), C_LOC(target_q), &
             C_LOC(source_x), C_LOC(source_y), C_LOC(source_z), C_LOC(source_q), &
             C_LOC(source_w), C_LOC(potential), &
             RBS_V, num_kernel_params, C_LOC(kernel_params), &
             singularity, approximation, compute_type, theta, degree, &
             max_source_leaf, max_target_leaf, size_check, beta, &
             verbosity)

    PRINT *, "RBS v total potential is: ", SUM(potential)

    DEALLOCATE(target_x, target_y, target_z, target_q, potential, source_x, source_y, &
             source_z, source_q, source_w);

    CALL MPI_FINALIZE(ierr)

END PROGRAM
