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

   A work-in-progress library for fast computation of N-body interactions on multiple GPUs,
   BaryTree implements barycentric Lagrange and Hermite polynomial interpolation treecodes.
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
on building and installing, see __INSTALL.md__ in this directory.


Examples
--------
See the __examples__ directory for several example executables that use __BaryTree__
and the Trilinos __Zoltan__ library for load balancing, in addition to an example 
using the library's interface for C programs. See __examples/README.md__ for more
details.


Interfaces
----------
See the __interfaces__ directory for __BaryTree__ interfaces for non-C programs.
Currently, there is a Python interface and an example script using that interface.
See __interfaces/README.md__ for more details.


References
----------
   Please refer to the following references for more background:
        
   - N. Vaughn, L. Wilson, and R. Krasny, A GPU-accelerated barycentric 
            Lagrange treecode, submitted to _Proc. 21st IEEE Int.
	    Workshop Parallel Distrib. Sci. Eng. Comput._ (PDSEC 2020) 
	    (2020).
	    
   - R. Krasny and L. Wang, A treecode based on barycentric Hermite 
            interpolation for electrostatic particle interactions,
	    _Comput. Math. Biophys._ __7__ (2019), 73-84.
		
   - H. A. Boateng and R. Krasny, Comparison of treecodes for
            computing electrostatic potentials in charged particle 
	    systems with disjoint targets and sources,
            _J. Comput. Chem._ __34__ (2013), 2159-2167.	
	   
   - J.-P. Berrut and L. N. Trefethen, Barycentric Lagrange interpolation,
            _SIAM Rev._ __46__ (2004), 501-517.

   - Z.-H. Duan and R. Krasny, An adaptive treecode for computing
            nonbonded potential energy in classical molecular systems,
            _J. Comput. Chem._ __22__ (2001), 184–195.

                                                    
License
-------
Copyright © 2019-2020, The Regents of the University of Michigan. Released under the [MIT License](LICENSE).
