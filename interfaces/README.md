Interfaces
----------

This folder contains interfaces between __BaryTree__ and other languages.
The treecode library itself contains `treedriverWrapper`. This function takes as input 
pointers to the particle arrays, as well runtime parameters such as kernel information, 
MAC parameter, and batch and cluster size. This function first constructs the particle 
structs then calls the `treedriver`.

The interfaces contained in these subdirectories are responsible for supplying 
`treedriverWrapper` with the necessary pointers to particle arrays and the runtime metadata.

----------

### Python

The python folder contains __treecodeWrappers.py__, which uses the ctypes module to load
the library, set the argument types, construct pointers to the numpy arrays, and call the
`treedriverWrapper`. 

__testTreecodeWrapper.py__ imports the python wrapper, generates 
some random particles, and calls the treecode once.

The `w` array is for quadrature weights when computing discrete convolution sums; 
it is set to ones for particle simulations.
