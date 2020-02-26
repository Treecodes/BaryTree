Examples
--------

This examples folder builds six executables:

1. `random_cube_cpu` and `random_cube_gpu`
2. `random_cube_large_cpu` and `random_cube_large_gpu`
3. `test_C_wrapper_cpu` and `test_C_wrapper_gpu`

- - -

### `random_cube` and `random_cube_large`


All of the random cube examples demonstrate the treecode's performance
using uniformly distributed cubes of random particles, load balanced
with Zoltan's recursive coordinate bisection.

The argument given to the executable is a parameter file that
specifies the run. An example is given here as `example.in`. For
example, one would run:

    random_cube_cpu example.in

to run the `random_cube_cpu` example with the parameters specified in
the file `example.in`.

The parameters in the infile are as follows:



- The `random_cube` examples are designed for reproducibility
of results–given a total number of particles across all ranks, the
actual random particles will be the same no matter how many ranks
are used.

- The `random_cube_large` examples are designed to test the
problem size limits of the treecode by overcoming limits in Zoltan's
maximum array sizes. Unlike the `random_cube` examples, which first 
generate all random particles and then use Zoltan to load balance them,
these examples generate a small number of particles, load balances
them, determines the resulting bounding boxes, and then generates the
specified number of random particles in those bounding boxes. The results
produced in terms of performanc and accuracy should be very similar to
the `random_cube` examples.

- - -

### `test_C_wrapper`

The `test_C_wrapper` examples demonstrate how to use the C
wrapper for the treecode.
