     ____                _______            
    |  _ \              |__   __|           
    | |_) | __ _ _ __ _   _| |_ __ ___  ___ 
    |  _ < / _` | '__| | | | | '__/ _ \/ _ \
    | |_) | (_| | |  | |_| | | | |  __|  __/
    |____/ \__,_|_|   \__, |_|_|  \___|\___|
                       __/ |                
                      |___/         
BaryTree
========

   This file is the partial documentation for BaryTree, a work-in-progress
   library for fast computation of N-body interactions on multiple GPUs.
   The current code employs an OpenACC GPU implementation.


   Authors:  
   - Leighton W. Wilson  (lwwilson@umich.edu) 
   - Nathan J. Vaughn  (njvaughn@umich.edu) 
   
   Department of Mathematics,
   University of Michigan, Ann Arbor.
   


Building
--------
This project uses CMake to manage and configure its build system. In principle, 
building this project is as simple as executing the following from the top level
directory of BaryTree:

    mkdir build; cd build; export CC=<C compiler>; cmake ..; make

Compiling GPU versions requires that a PGI C compiler be used. For more information
on building and installing, see the INSTALL.md file in this directory.


References
----------
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


                                                     
License
-------
Copyright © 2019-2020, The Regents of the University of Michigan. Released under the [MIT License](LICENSE).
