#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "../utilities/timers.h"
#include "../particles/struct_particles.h"
#include "../run_params/struct_run_params.h"
#include "../interaction_compute/interaction_compute.h"

#include "directdriver.h"


void directdriver(struct Particles *sources, struct Particles *targets, struct RunParams *run_params,
                  double *potential, double *time_direct)
{
    time_direct[0] = 0.0;
    time_direct[1] = 0.0;
    time_direct[2] = 0.0;
    time_direct[3] = 0.0;

    START_TIMER(&time_direct[2]);
    InteractionCompute_Direct(potential, sources, targets, run_params);
    STOP_TIMER(&time_direct[2]);

    return;

} /* END of function directdriver */
