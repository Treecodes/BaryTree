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
   

   Please include the following references in any work that utilizes this code:
		
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
  tree.exe:     Runs the particle-cluster or cluster-particle treecode for Coulomb or
                screened Coulomb interactions.
	     
  direct.exe:   Directly computes Coulomb or screened Coulomb interactions.
  
  txt2bin.exe:  Turns a text file of sources (pqr or x/y/z/q file) or targets into a 
  		binary file for use by the tree.exe or direct.exe executables.


                     
Building
------------------------------
To compile with GPU support, first load all modules necessary for PGI compilers.
For instance, on Michigan Flux, first run:
	      
	module load cuda cupti pgi openmpi/1.10.2/pgi/16.4
	      
To profile the code, either add the flag `-ta=nvidia,time` to the Makefile in
src, or run the following commands:
	      
	export PGI_ACC_TIME=1
	export ACC_NOTIFY=1
	      
The code is compiled with a `make` command in the src directory. This will generate
executables in the bin directory. By default, on Linux, the compilation will rely
on PGI compilers for GPU support, and by default, on Mac, the compilation will not
include GPU support. Edit the flags in Makefile to change this behavior.
   
   
              
                                                     
Input for tree.exe:
-------------------------------

      Enter the following as command line arguments:
 
                   sampin1:  sources input file 
                   sampin2:  targets input file 
                   sampin3:  direct calc potential input file 
                   sampout:  tree calc potential output file
                  numparsS:  number of sources 
                  numparsT:  number of targets
                     theta:  multipole acceptance criterion
                     order:  order of treecode Taylor expansion 
                 tree_type:  0--cluster-particle, 1--particle-cluster 
                maxparnode:  maximum particles in leaf 
                     kappa:  screened Coulomb parameter
                  pot_type:  0--Coulomb, 1--screened Coulomb
                     pflag:  distribute 0--targets, 1--sources
                     sflag:  on distr 0--sort, 1--no sort
                     dflag:  if sorted, direction 0--x, 1--y, 2--z	
		     batch:  maximum particles in batch

License
-------
Copyright © 2019, The Regents of the University of Michigan. Released under the [MIT License](LICENSE).
