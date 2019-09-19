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

                     
License
-------
Copyright © 2019, The Regents of the University of Michigan. Released under the [MIT License](LICENSE).
