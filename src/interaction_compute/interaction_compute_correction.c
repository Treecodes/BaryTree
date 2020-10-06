#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "../particles/struct_particles.h"
#include "../run_params/struct_run_params.h"

#include "../kernels/coulomb/coulomb.h"
#include "../kernels/yukawa/yukawa.h"
#include "../kernels/regularized-coulomb/regularized-coulomb.h"
#include "../kernels/regularized-yukawa/regularized-yukawa.h"
#include "../kernels/atan/atan.h"
#include "../kernels/sin-over-r/sin-over-r.h"
#include "../kernels/mq/mq.h"
#include "../kernels/user_kernel/user_kernel.h"

#include "interaction_compute.h"


void InteractionCompute_SubtractionPotentialCorrection(double *potential,
                           struct Particles *targets, struct RunParams *run_params)
{
    int num_targets  = targets->num;
    double *target_q = targets->q;

    if (run_params->singularity == SUBTRACTION) {
        if (run_params->kernel == COULOMB) {
            K_Coulomb_SS_Correction(potential, target_q, num_targets, run_params);

        } else if (run_params->kernel == REGULARIZED_COULOMB) {
            K_RegularizedCoulomb_SS_Correction(potential, target_q, num_targets, run_params);

        } else if (run_params->kernel == YUKAWA) {
            K_Yukawa_SS_Correction(potential, target_q, num_targets, run_params);

        } else if (run_params->kernel == REGULARIZED_YUKAWA) {
            K_RegularizedYukawa_SS_Correction(potential, target_q, num_targets, run_params);

        }
    }

    return;
}
