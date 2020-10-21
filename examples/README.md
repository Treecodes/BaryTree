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
| `distribution`    | Underlying particle distribution: `UNIFORM`, `GAUSSIAN`, `EXPONENTIAL`, `PLUMMER`, or `PLUMMER_SYMMETRIC`.
| `degree`          | Degree of polynomial interpolation. 
| `theta`           | Multipole acceptance criterion (MAC).
| `max_per_source_leaf` | Maximum number of particles per source tree leaf (or source batch, for `CLUSTER_PARTICLE`).
| `max_per_target_leaf` | Maximum number of particles per target tree leaf (or target batch, for `PARTICLE_CLUSTER`).
| `beta`            | Automatic tuning accuracy parameter. Number in [0,1], higher is more accurate. 
| `compute_type`    | Type of treecode method. `CLUSTER_PARTICLE`, `PARTICLE_CLUSTER` (i.e. BLTC), `CLUSTER_CLUSTER` (i.e. BLDTT).
| `approximation`   | Type of polynomial: `LAGRANGE` and `HERMITE`. `HERMITE` is incompatible with cluster-cluster.
| `kernel_name`     | Name of interaction kernel: `COULOMB`, `YUKAWA`, `REGULARIZED_COULOMB`, `REGULARIZED_YUKAWA`, `SIN_OVER_R`, `USER`.
| `kernel_params`   | Comma separated list of parameters for given kernel.
| `run_direct`      | Run direct calculation for error comparison: `ON` or `OFF`.
| `verbosity`       | Determines verbosity level of output. Integer `0`, `1`, `2`, `3`. Higher means more output.
| `slice`           | Determines the proportion of target sites at which the direct calculation is performed for error comparison. 10 would mean every 10th target is sampled.


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
