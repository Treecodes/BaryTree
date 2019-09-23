Examples
========

Making example point sets
-------------------------
Use the make\_points.py script to generate example source and target binary files for running with the executables.
For example, to generate a file of sources with $N$ particles, run:

    python3 make_points.py sources sources_example.bin N

To generate a file of targets with $M$ targets, run:

    python3 make_points.py targets targets_example.bin M

Using the library
-----------------
`example.c` demonstrates using the library in another program. This program can be built by configuring CMake
with the flag `-DBUILD\_EXAMPLES=ON` when building the executables and libraries in the top level directory.

Using the executables
---------------------
Two short bash scripts, `exec_example_cpu.sh` and `exec_example_gpu.sh`, are included to demonstrate generating
 example point sets, running the direct executable to produce
a benchmark solution, and then running the tree executable to generate an approximate solution. What follows is a
description of what the `exec_example_cpu.sh` script does: 

First produce example point sets using the make_points.py script:

    python3 make_points.py sources sources_1E4.bin 10000
    python3 make_points.py targets targets_1E4.bin 10000

Assuming the executables have been installed, run the direct executable to generate a reference solution for the Coulomb
potential, using 1 OpenMP thread:

    direct-cpu sources_1E4.bin targets_1E4.bin direct_benchmark.bin direct_summary.csv 10000 10000 0.0 0 1

Now run the tree executable using Lagrange barycentric interpolation for the Coulomb potential, with a MAC of 0.5 and
an interpolation order of 5, using 1 OpenMP thread, and compare it to the reference solution:

    tree-cpu sources_1E4.bin targets_1E4.bin direct_benchmark.bin tree_summary.csv 10000 10000 0.5 5 500 5 0 0.0 1


                     
License
-------
Copyright © 2019, The Regents of the University of Michigan. Released under the [MIT License](LICENSE).
