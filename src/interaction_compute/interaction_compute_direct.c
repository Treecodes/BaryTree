#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "../utilities/array.h"

#include "../particles/struct_particles.h"
#include "../run_params/struct_run_params.h"

#include "../kernels/coulomb/coulomb.h"
#include "../kernels/yukawa/yukawa.h"
#include "../kernels/regularized-coulomb/regularized-coulomb.h"
#include "../kernels/regularized-yukawa/regularized-yukawa.h"
#include "../kernels/atan/atan.h"
#include "../kernels/sin-over-r/sin-over-r.h"
#include "../kernels/mq/mq.h"

#include "interaction_compute.h"


void InteractionCompute_Direct(double *potential,
                               struct Particles *sources, struct Particles *targets,
                               struct RunParams *run_params)
{

    int num_sources   = sources->num;
    double *source_x  = sources->x;
    double *source_y  = sources->y;
    double *source_z  = sources->z;
    double *source_q  = sources->q;
    double *source_w  = sources->w;

    int num_targets   = targets->num;
    double *target_x  = targets->x;
    double *target_y  = targets->y;
    double *target_z  = targets->z;
    double *target_q  = targets->q;

#ifdef OPENACC_ENABLED
    #pragma acc data copyin(source_x[0:num_sources], source_y[0:num_sources], source_z[0:num_sources], \
                            source_q[0:num_sources], source_w[0:num_sources], \
                            target_x[0:num_targets], target_y[0:num_targets], target_z[0:num_targets], \
                            target_q[0:num_targets]), copy(potential[0:num_targets])
#endif
    {

/**********************************************************/
/**************** COMPLETE DIRECT SUM *********************/
/**********************************************************/


    /***************************************/
    /********* Coulomb *********************/
    /***************************************/

    if (run_params->kernel == COULOMB) {

        if (run_params->singularity == SKIPPING) {

            K_Coulomb_Direct(num_targets, num_sources, 0, 0,
                    target_x, target_y, target_z,
                    source_x, source_y, source_z, source_q, source_w,
                    run_params, potential, 0);

        } else if (run_params->singularity == SUBTRACTION) {

            K_Coulomb_SS_Direct(num_targets, num_sources, 0, 0,
                    target_x, target_y, target_z, target_q,
                    source_x, source_y, source_z, source_q, source_w,
                    run_params, potential, 0);

        } else {
            printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
            exit(1);
        }

    /***************************************/
    /********* Yukawa **********************/
    /***************************************/

    } else if (run_params->kernel == YUKAWA) {

        if (run_params->singularity == SKIPPING) {

            K_Yukawa_Direct(num_targets, num_sources, 0, 0,
                    target_x, target_y, target_z,
                    source_x, source_y, source_z, source_q, source_w,
                    run_params, potential, 0);

        } else if (run_params->singularity == SUBTRACTION) {

            K_Yukawa_SS_Direct(num_targets, num_sources, 0, 0,
                    target_x, target_y, target_z, target_q,
                    source_x, source_y, source_z, source_q, source_w,
                    run_params, potential, 0);

        } else {
            printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
            exit(1);
        }

    /***************************************/
    /********* Regularized-Coulomb *********/
    /***************************************/

    } else if (run_params->kernel == REGULARIZED_COULOMB) {

        if (run_params->singularity == SKIPPING) {

            K_RegularizedCoulomb_Direct(num_targets, num_sources, 0, 0,
                    target_x, target_y, target_z,
                    source_x, source_y, source_z, source_q, source_w,
                    run_params, potential, 0);

        } else if (run_params->singularity == SUBTRACTION) {

            K_RegularizedCoulomb_SS_Direct(num_targets, num_sources, 0, 0,
                    target_x, target_y, target_z, target_q,
                    source_x, source_y, source_z, source_q, source_w,
                    run_params, potential, 0);

        }


    /***************************************/
    /********* Regularized-Yukawa **********/
    /***************************************/

    } else if (run_params->kernel == REGULARIZED_YUKAWA) {

        if (run_params->singularity == SKIPPING) {

            K_RegularizedYukawa_Direct(num_targets, num_sources, 0, 0,
                        target_x, target_y, target_z,
                        source_x, source_y, source_z, source_q, source_w,
                        run_params, potential, 0);

        } else if (run_params->singularity == SUBTRACTION) {

            K_RegularizedYukawa_SS_Direct(num_targets, num_sources, 0, 0,
                        target_x, target_y, target_z, target_q,
                        source_x, source_y, source_z, source_q, source_w,
                        run_params, potential, 0);
        }


    /***************************************/
    /********* Atan ************************/
    /***************************************/

    } else if (run_params->kernel == ATAN) {

            K_Atan_Direct(num_targets, num_sources, 0, 0,
                        target_x, target_y, target_z,
                        source_x, source_y, source_z, source_q, source_w,
                        run_params, potential, 0);


    /***************************************/
    /********* Sin Over R ******************/
    /***************************************/

    } else if (run_params->kernel == SIN_OVER_R) {

            K_SinOverR_Direct(num_targets, num_sources, 0, 0,
                        target_x, target_y, target_z,
                        source_x, source_y, source_z, source_q, source_w,
                        run_params, potential, 0);

    /***************************************/
    /************ MQ ***********************/
    /***************************************/

    } else if (run_params->kernel == MQ) {

            K_MQ_Direct(num_targets, num_sources, 0, 0,
                        target_x, target_y, target_z,
                        source_x, source_y, source_z, source_q, source_w,
                        run_params, potential, 0);

    } else {
        printf("**ERROR** INVALID KERNEL. EXITING.\n");
        exit(1);
    }

#ifdef OPENACC_ENABLED
        #pragma acc wait
#endif
    } // end acc data region

    return;
}
