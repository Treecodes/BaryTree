#!/bin/bash

echo BEGINNING CPU TESTS
tests_cpu
mpirun -np 2 tests_cpu_mpi

echo BEGINNING GPU TESTS
tests_gpu