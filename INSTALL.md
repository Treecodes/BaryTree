Installing
----------

In a best case scenario, building and installing the libraries and examples should be as simple as this:

    mkdir build; cd build; export CC=<C compiler>;
    cmake .. -DCMAKE_INSTALL_PREFIX=<install location>;
    make -j install;

This assumes that you have a few things:
1. a sane C compiler,
2. a sane MPI installation that agrees with your C compiler,
3. CMake version 3.9 or newer,
4. an installed Trilinos Zoltan library (for building examples).

If you don't have an install of Zoltan, then you must turn off building of
examples with `-DBUILD_EXAMPLES=OFF`.

Compiling GPU versions requires that a PGI C compiler be used. If another compiler
other than pgcc is used, for instance gcc or icc, support for building GPU versions
will be automatically turned off during configuration.

CMake Flags
-----------
A list of the most useful CMake flags to use during configure:

 - `-DCMAKE_RELEASE_TYPE={Debug, Release}`      build either the debug or release version
 - `-DENABLE_GPU_BUILD={ON, OFF}`               toggle whether to build the GPU versions
 - `-DCMAKE_INSTALL_PREFIX=<where to install>`  specify install location for `make install`
 - `-DBUILD_EXAMPLES={ON, OFF}`                 toggle whether to build examples (requires Zoltan)
 - `-DBUILD_SHARED_LIBS={ON, OFF}`              toggle whether to build libraries as shared or static objects
 - `-DZoltan_DIR=<location of Zoltan cmake>`    specify location of Zoltan CMake configuration file if not picked up by CMake automatically (typically `lib/cmake/Zoltan` of wherever Trilinos was installed)
 
 If the Zoltan install isn't picked up automatically, you can also add the install location of Trilinos or Zoltan to the CMake module search path with `-DCMAKE_PREFIX_PATH=<location of Zoltan install>`. This is an alternative to explicitly setting `Zoltan_DIR`.
    
Testing
-------
After building, you can test the build by running `ctest` or `make test` from the build
directory. This performs a series of simple serial tests.
