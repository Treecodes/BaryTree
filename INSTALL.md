Installing
----------

In a best case scenario, building and installing the libraries and examples should be as simple as this:

    export CC=<YOUR C COMPILER HERE>
    mkdir build;
    cd build;
    cmake .. -DCMAKE_INSTALL_PREFIX=<YOUR INSTALL LOCATION HERE>
    make -j install

This assumes a few things. One, that you have a sane C compiler; two, that you have a sane MPI installation
that agrees with your C compiler; three, that you have an installed Zoltan library (for building examples and tests);
and four, that you have a CMake version 3.9 or newer. CMake should be able to detect where your MPI and Zoltan
libraries are located.

If you don't have an install of Zoltan, and don't want to build the examples and tests, then add the following
flags to the `cmake` command: -DBUILD_EXAMPLES=OFF -DBUILD_TESTS=OFF

After building, you can test the build by running `ctest` from the build directory.

Note that building with GPU support requires PGI compilers. If the compilers are not PGI, then GPU support will be
automatically turned off.
