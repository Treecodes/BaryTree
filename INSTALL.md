Installing
----------

In a best case scenario, building and installing the libraries and examples should be as simple as this:

    mkdir build; cd build; export CC=<C compiler>;
    cmake .. -DCMAKE_INSTALL_PREFIX=<install location>;
    make -j install;

This assumes that you have a few things:
1. a sane C compiler,
2. a sane MPI installation that agrees with your C compiler,
3. an installed Zoltan library (for building examples and tests),
4. CMake version 3.9 or newer. CMake should be able to detect where your MPI and Zoltan libraries are located.


If you don't have an install of Zoltan, and don't want to build the examples and tests, then add the following
flags to the `cmake` command: `-DBUILD_EXAMPLES=OFF -DBUILD_TESTS=OFF`


After building, you can test the build by running `ctest` from the build directory.


Note that building with GPU support requires PGI compilers. If the compilers are not PGI, then GPU support will be
automatically turned off.
