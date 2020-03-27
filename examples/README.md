Examples
--------

This examples folder builds six executables:

1. __random_cube_cpu__ and __random_cube_gpu__
2. __random_cube_reproducible_cpu__ and __random_cube_reproducible_gpu__
3. __testBaryTreeInterface_cpu__ and __testBaryTreeInterface_gpu__

- - -

#### __random_cube__ and __random_cube_reproducible__

All of the random cube examples demonstrate the treecode's performance
using a cube of uniformly distributed random particles, load balanced
with Zoltan's recursive coordinate bisection.

The argument given to the executable is a parameter file that
specifies the run. An example is given here as __example.in__. For
example, one would run:

    mpirun -n 2 random_cube_cpu example.in

to run the __random_cube_cpu__ example with the parameters specified in
the file __example.in__ across two ranks.

The parameters that can be specified in the infile are as follows:
| Parameter         | Description
|-------------------|------------------
| `num_particles`   | Number of sources and targets. Its use is exclusive with the `num_sources` and `num_targets` parameters.
| `num_sources`     | Number of sources.
| `num_targets`     | Number of targets.
| `order`           | Order of polynomial interpolation. 
| `theta`           | Multipole acceptance criterion (MAC).
| `max_per_leaf`    | Maximum number of particles per tree leaf.
| `max_per_batch`   | Maximum number of particles per batch.
| `kernel_name`     | Name of interaction kernel: `yukawa` or `coulomb`.
| `approximation`   | Type of polynomial: `lagrange` and `hermite`. 
| `size_check`      | If the product of this parameter and the number of interpolation points in a cluster is greater than the number of particles in the cluster, then the interaction will be performed directly even if the MAC is accepted.
| `run_direct`      | Run direct calculation for error comparison: `on` or `off`
| `verbosity`       | Determines verbosity level of output. `0` is quiet, `1` is verbose.
| `slice`           | Determines the proportion of target sites at which the direct calculation is performed for error comparison.
| `kernel_params`   | Comma separated list of parameters for given kernel.

Note the difference between these executables:

- The __random_cube__ examples are designed to test the
problem size limits of the treecode by overcoming limits in Zoltan's
maximum array sizes. Unlike the __random_cube_reproducible__ examples, which first 
generate all random particles and then use Zoltan to load balance them,
these examples generate a small number of particles, load balances
them, determines the resulting bounding boxes, and then generates the
specified number of random particles in those bounding boxes. The results
produced in terms of performance and accuracy should be very similar to
the __random_cube_reproducible__ examples.

- The __random_cube_reproducible__ examples are designed for reproducibility
of results. Given a total number of particles across all ranks, the
actual random particles will be the same no matter how many ranks
are used (given that the executable is run on the same computational
resource). Additionally, this example requires that the number of sources
and targets be equal.

- - -

#### __testBaryTreeInterface__

The __testBaryTreeInterface__ examples demonstrate how to use the C wrapper 
for the treecode. A C program that links to the __BaryTree__ library can, 
in fact, directly use the `treedriver` function if the calling program 
implements the particle and kernel struct used by `treedriver` 
(as done in the above examples). The `BaryTreeInterface` function, 
however, takes source and target particle arrays directly.
