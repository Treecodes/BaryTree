BaryTree
========

   This file is the partial documentation for hybrid-gpu-treecode, 
   a work-in-progress set of routines for approximating the electrostatic 
   energy at a set of M targets due to a disjoint set of N source particles.
   The current code employs an OpenACC GPU implementation to handle direct
   source-target interactions, and can also be compiled with MPI 
   capabilities.


   Authors:  
   - Leighton W. Wilson  (lwwilson@umich.edu) 
   - Nathan J. Vaughn  (njvaughn@umich.edu) 
   
   Department of Mathematics,
   University of Michigan, Ann Arbor
   

   Please refer to the following references for more background:
		
   - Boateng. H. A., Krasny, R.: Comparison of Treecodes for
            Computing Electrostatic Potentials in Charged Particle 
	    Systems with Disjoint Targets and Sources.
            _J. Comput. Chem._ (2013)	 

   - Duan, Z.-H., Krasny, R.: An adaptive treecode for computing
            nonbonded potential energy in classical molecular systems.
            _J. Comput. Chem._ __22__ (2001) 184–195
 
   - Lindsay, K., Krasny, R.: A particle method and adaptive treecode
            for vortex sheet motion in 3-D flow. _J. Comput. Phys._ __172__
            (2001) 879–907

   - Deng, Q., Driscoll, T. A: A Fast Treecode for Multiquadric 
            Interpolation with Varying Shape Parameters.
            _SIAM J. Sci. Comput._ __34__ (2012) A1126–A1140



Summary of executables
----------------------
- tree-cpu:    Runs the particle-cluster Lagrange or Hermite barycentric
              treecode for Coulomb or screened Coulomb interactions on
              CPUs, using OpenMP for single-node parallelization.
	     
- tree-gpu:    Runs the particle-cluster Lagrange or Hermite barycentric
              treecode for Coulomb or screened Coulomb interactions on
              one or more GPUs connected to one compute node.
	     
- direct-cpu:  Directly computes Coulomb or screened Coulomb interactions
              on CPUs, using OpenMP for single-node parallelization.

- direct-gpu:  Directly computes Coulomb or screened Coulomb interactions
              on one or more GPUs connected to one compute node.
  

                     
Building
------------------------------
This project uses CMake to manage and configure its build system. In principle, 
building this project is as simple as executing the following from the top level
directory of BaryTree:

    mkdir build; cd build; export CC=<C compiler>; cmake ..; make

Compiling GPU versions requires that a PGI C compiler be used. If another compiler
other than pgcc is used, for instance gcc or icc, support for building GPU versions
will be automatically turned off during configuration.

Some potentiall useful CMake flags during configure:

    -DCMAKE_RELEASE_TYPE={Debug, Release}   build either the debug or release version
    -DENABLE_GPU_BUILD={ON, OFF}   manually toggle whether to build the GPU versions
    -DCMAKE_INSTALL_PREFIX=/where/to/install   specify install location for `make install`
    

    
	      
   
   
              
                                                     
Running the executables
-----------------------
To run `direct-cpu` or `direct-gpu`, enter the following as command line arguments:

              infile 1:  sources input file 
              infile 2:  targets input file 
              infile 3:  direct calc potential binary output file 
            csv output:  results summary to CSV file
              numparsS:  number of sources 
              numparsT:  number of targets 
                 kappa:  screened Coulomb parameter 
              pot type:  0--Coulomb
                         1--screened Coulomb 
      num threads/GPUs:  number of OpenMP threads or available GPUs

Running `direct-cpu --help` or `direct-gpu --help` help will produce a list of these command
line arguments

To run `tree-cpu` or `tree-gpu`, enter the following as command line arguments:

              infile 1:  sources binary input file
              infile 2:  targets binary input file
              infile 3:  direct calc potential binary input file 
            csv output:  results summary to CSV file 
              numparsS:  number of sources 
              numparsT:  number of targets 
                 theta:  multipole acceptance criterion 
                 order:  number of Chebyshev interp. pts per Cartesian direction 
             leaf size:  maximum particles in leaf 
            batch size:  maximum size of target batch 
       pot/approx type:  0--Coulomb, Lagrange approx.
                         1--screened Coulomb/Yukawa, Lagrange approx.
                         4--Coulomb, Hermite approx.
                         5--screened Coulomb/Yukawa, Hermite approx.
                 kappa:  screened Coulomb parameter 
      num threads/GPUs:  number of OpenMP threads or available GPUs

For example, running:

    tree-gpu sources.bin targets.bin direct_result.bin tree_result.csv 100000 100000 0.5 5 500 500 4 0 2

would run tree-gpu on 2 GPUs for the Coulomb potential, using Hermite interpolation with
5 interpolation points per Cartesian direction and a MAC of 0.5. The executable would 
compare the results to a direct result produced by direct-gpu or direct-cpu and saved in
the file direct\_result.bin. 

Running `tree-cpu --help` or `tree-gpu --help` will produce a list of these command line arguments.
Running with a non-existent file for the direct calc potential binary input file argument (the third
argument) will still produce timing info, just with no benchmark comparison.



License
-------
Copyright © 2019, The Regents of the University of Michigan. Released under the [MIT License](LICENSE).
