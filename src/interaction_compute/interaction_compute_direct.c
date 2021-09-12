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
#include "../kernels/tcf/tcf.h"
#include "../kernels/dcf/dcf.h"

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

    double target_xdd = targets->xdd;
    double target_ydd = targets->ydd;
    double target_zdd = targets->zdd;

    double target_xmin = targets->xmin;
    double target_ymin = targets->ymin;
    double target_zmin = targets->zmin;

    double target_xmax = targets->xmax;
    double target_ymax = targets->ymax;
    double target_zmax = targets->zmax;

    int target_xdim = targets->xdim;
    int target_ydim = targets->ydim;
    int target_zdim = targets->zdim;


#ifdef OPENACC_ENABLED
//    #pragma acc data copyin(source_x[0:num_sources], source_y[0:num_sources], source_z[0:num_sources], \
//                            source_q[0:num_sources], source_w[0:num_sources], \
//                            target_x[0:num_targets], target_y[0:num_targets], target_z[0:num_targets], \
//                            target_q[0:num_targets]), copy(potential[0:num_targets])
    #pragma acc data copyin(source_x[0:num_sources], source_y[0:num_sources], source_z[0:num_sources], \
                            source_q[0:num_sources], source_w[0:num_sources]), copy(potential[0:num_targets])
#endif
    {

/* * ********************************************************/
/* * ************** COMPLETE DIRECT SUM *********************/
/* * ********************************************************/


    /* * *************************************/
    /* * ******* Coulomb *********************/
    /* * *************************************/

    if (run_params->kernel == COULOMB) {

        if (run_params->singularity == SKIPPING) {

            K_Coulomb_Direct(num_targets, num_sources, 0, 0,
                    target_x, target_y, target_z,
                    source_x, source_y, source_z, source_q, source_w,
                    run_params, potential, 0);

        } else {
            printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
            exit(1);
        }

    /* * *************************************/
    /* * ******* Yukawa **********************/
    /* * *************************************/

    } else if (run_params->kernel == YUKAWA) {

        if (run_params->singularity == SKIPPING) {

            K_Yukawa_Direct(num_targets, num_sources, 0, 0,
                    target_x, target_y, target_z,
                    source_x, source_y, source_z, source_q, source_w,
                    run_params, potential, 0);

        } else {
            printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
            exit(1);
        }


    /* * *************************************/
    /* * ******* TCF *************************/
    /* * *************************************/

    } else if (run_params->kernel == TCF) {
            K_TCF_Direct(0,  target_xdim-1,
                         0,  target_ydim-1,
                         0,  target_zdim-1,
                         target_xmin, target_ymin, target_zmin,
                         
                         target_xdd,  target_ydd,  target_zdd,
                         target_xdim, target_ydim, target_zdim,

                         num_sources, 0,
                         source_x, source_y, source_z, source_q,

                         run_params, potential, 0);

            //K_TCF_Direct(num_targets, num_sources, 0, 0,
            //            target_x, target_y, target_z,
            //            source_x, source_y, source_z, source_q, source_w,
            //            run_params, potential, 0);
                        
                        
    /* * *************************************/
    /* * ******* DCF *************************/
    /* * *************************************/

    } else if (run_params->kernel == DCF) {

            K_DCF_Direct(num_targets, num_sources, 0, 0,
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
