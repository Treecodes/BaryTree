#!/bin/bash



mpirun -np 1 zoltan_example_gpu  500000 7 0.8 2000 2000 coulomb 0.0 skipping lagrange 1 0 0
mpirun -np 1 zoltan_example_gpu 1000000 7 0.8 2000 2000 coulomb 0.0 skipping lagrange 1 0 0
mpirun -np 1 zoltan_example_gpu 4000000 7 0.8 2000 2000 coulomb 0.0 skipping lagrange 1 0 0 