#include <math.h>
#include <float.h>
#include <stdio.h>

#include "../../struct_kernel.h"
#include "regularized-coulomb_ss_correction.h"


void K_RegularizedCoulomb_SS_Correction(double *potential, double *target_q,
                                        int numTargets, struct kernel *kernel)
{
    double alpha = kernel->parameters[0];
    double param = 2.0 * M_PI * alpha * alpha;
    for (int i = 0; i < numTargets; i++) potential[i] += param * target_q[i];

    return;
}
