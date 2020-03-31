Adding New Kernels
------------------

Steps for adding a new kernel named `custom-kernel` with support for particle-cluster are described below. Adding support for cluster-cluster and cluster-particle follows similarly. Consult existing kernel directories for more details.

1. Add the name of the new kernel to the end of the KERNEL enum in `src/utilities/enums.h`. If you plan to use the Python wrapper, add it to the Kernel class in `interfaces/python/BaryTreeInterface.py` as well.

2. Create a new directory `custom-kernel` in `src/kernels/`.

3. Create new source files (`custom-kernel_direct.c`, `custom-kernel_direct.h`, `custom-kernel_pc.c`, `custom-kernel_pc.h`) in `src/kernels/custom-kernel` containing three batch-cluster interaction functions: 
	- `K_CustomKernel_Direct( )`
	- `K_CustomKernel_PC_Lagrange( )`
	- `K_CustomKernel_PC_Lagrange( )`
	
4. Create a new source file `custom-kernel.h` in `src/kernels/custom-kernel/`. `#include` in this file all other headers associated with this kernel (`custom-kernel_direct.h` and `custom-kernel_pc.h`).

5. Edit `interaction_compute_pc.c`:
	1. Include `custom-kernel.h` in `interaction_compute_pc.c` (`#include "kernels/custom_kernel.h"`).
	2. Add your custom kernel in several places, following the format for the already-present kernels, like `if (run_params->kernel == KERNEL_NAME)`:
		- In the POTENTIAL FROM APPROX subsection, add your Lagrange and/or Hermite kernels
		- In the POTENTIAL FROM DIRECT subsection, add your direct interaction kernel
			
6. Edit `interaction_compute_direct.c`:
	1. Include `custom-kernel.h` in `interaction_compute_pc.c` (`#include "kernels/custom_kernel.h"`).
	2. Add your direct interaction kernel, following the format for the already-present kernels.	

7. Add your files to `src/CMakeLists.txt`:
	1. Add `custom-kernel_direct.c`, `custom-kernel_direct.h`, `custom-kernel_pc.c`, `custom-kernel_pc.h`, and `kernels/custom_kernel.h` to a new `SRCS_K_CUSTOM_KERNEL` list, following the format of the other kernels.
	2. Add `${SRCS_K_CUSTOM_KERNEL}` to the `SRCS_KERNELS` list.
