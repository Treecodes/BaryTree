Interfaces
----------

This folder contains interfaces between __BaryTree__ and other languages.
The BaryTree library itself contains `BaryTreeInterface`. This function takes as input 
pointers to the particle arrays, as well runtime parameters such as kernel information, 
MAC parameter, and batch and cluster size. This function first constructs the particle 
structs then calls the `treedriver`.

The interfaces contained in these subdirectories are responsible for supplying 
`BaryTreeInterface` with the necessary pointers to particle arrays and the runtime metadata.

----------

### Python

The Python folder contains __BaryTreeInterface.py__, which uses the `ctypes` module to load
the library, set the argument types, construct pointers to the `numpy` arrays, and call the
`BaryTreeInterface`. 

__testBaryTreeInterface.py__ imports the Python wrapper, generates some random particles, 
and calls the treecode once.

The `w` array is for quadrature weights when computing discrete convolution sums; 
it is set to ones for particle simulations.
