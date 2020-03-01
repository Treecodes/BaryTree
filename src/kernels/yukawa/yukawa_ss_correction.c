#include <math.h>
#include <float.h>
#include <stdio.h>

#include "../../struct_kernel.h"
#include "yukawa_ss_correction.h"


void K_Yukawa_SS_Correction(double *potential, double *target_q,
                            int numTargets, struct kernel *kernel)
{
    double kernel_parameter=kernel->parameters[0];
    double param = 4.0 * M_PI / kernel_parameter / kernel_parameter;
    for (int i = 0; i < numTargets; i++) potential[i] += param * target_q[i];

    return;
}
