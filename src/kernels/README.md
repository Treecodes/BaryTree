Adding New Kernels
------------------

Steps for adding a new kernel named `custom_kernel`:

1. Create new source files (__custom_kernel.c__, __custom_kernel.h__) in `src/kernels/` containing three batch-cluster interaction functions: 
	-`customKernelDirect( )`
	-`customKernelApproximationLagrange( )`
	-`customKernelApproximationHermite( )`

2. Edit __interaction_compute.c__
	1. Include __custom-kernel.h__ in __interaction_compute.c__ (i.e.  `#include "kernels/custom_kernel.h"`)
	2. Add your custom kernel in several places, following the format for the already-present kernels, using string comparison to distinguish your kernel, e.g. `if (strcmp(kernel->name, "custom_kernel") == 0)`
		- Inside of `Interaction_PC_Compute()`
			- In the POTENTIAL FROM APPROX subsection, add your Lagrange and/or Hermite kernels.
			- In the POTENTIAL FROM DIRECT subsection, add your direct interaction kernel.
		- Inside of `Interaction_Direct_Compute()`
			- Add direct interaction kernel for performing direct-sum reference calculations.

3. Add your files to `src/CMakeLists.txt`
	- add `kernels/custom-kernel.c` and `kernels/custom-kernel.h` to the `SRCS_TREEDRIVER` list.


