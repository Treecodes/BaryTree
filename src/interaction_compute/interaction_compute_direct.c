#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>

#include "../utilities/array.h"

#include "../tree/struct_nodes.h"
#include "../particles/struct_particles.h"
#include "../run_params/struct_run_params.h"

#include "../kernels/coulomb/coulomb.h"
#include "../kernels/yukawa/yukawa.h"
#include "../kernels/regularized-coulomb/regularized-coulomb.h"
#include "../kernels/regularized-yukawa/regularized-yukawa.h"
#include "../kernels/atan/atan.h"


#include "interaction_compute.h"


void InteractionCompute_Direct(double *source_x, double *source_y, double *source_z,
                               double *source_q, double *source_w,
                               double *target_x, double *target_y, double *target_z, double *target_q,
                               double *pointwisePotential, int numSources, int numTargets,
                               struct RunParams *run_params)
{


#ifdef OPENACC_ENABLED
    #pragma acc data copyin(source_x[0:numSources], source_y[0:numSources], source_z[0:numSources], \
                            source_q[0:numSources], source_w[0:numSources], \
                            target_x[0:numTargets], target_y[0:numTargets], target_z[0:numTargets], \
                            target_q[0:numTargets]), copy(pointwisePotential[0:numTargets])
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

            K_Coulomb_Direct(numTargets, numSources, 0, 0,
                    target_x, target_y, target_z,
                    source_x, source_y, source_z, source_q, source_w,
                    run_params, pointwisePotential, 0);

        } else if (run_params->singularity == SUBTRACTION) {

            K_Coulomb_SS_Direct(numTargets, numSources, 0, 0,
                    target_x, target_y, target_z, target_q,
                    source_x, source_y, source_z, source_q, source_w,
                    run_params, pointwisePotential, 0);

        } else {
            printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
            exit(1);
        }

    /***************************************/
    /********* Yukawa **********************/
    /***************************************/

    } else if (run_params->kernel == YUKAWA) {

        if (run_params->singularity == SKIPPING) {

            K_Yukawa_Direct(numTargets, numSources, 0, 0,
                    target_x, target_y, target_z,
                    source_x, source_y, source_z, source_q, source_w,
                    run_params, pointwisePotential, 0);

        } else if (run_params->singularity == SUBTRACTION) {

            K_Yukawa_SS_Direct(numTargets, numSources, 0, 0,
                    target_x, target_y, target_z, target_q,
                    source_x, source_y, source_z, source_q, source_w,
                    run_params, pointwisePotential, 0);

        } else {
            printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
            exit(1);
        }

    /***************************************/
    /********* Regularized-Coulomb *********/
    /***************************************/

    } else if (run_params->kernel == REGULARIZED_COULOMB) {

        if (run_params->singularity == SKIPPING) {

            K_RegularizedCoulomb_Direct(numTargets, numSources, 0, 0,
                    target_x, target_y, target_z,
                    source_x, source_y, source_z, source_q, source_w,
                    run_params, pointwisePotential, 0);

        } else if (run_params->singularity == SUBTRACTION) {

            K_RegularizedCoulomb_SS_Direct(numTargets, numSources, 0, 0,
                    target_x, target_y, target_z, target_q,
                    source_x, source_y, source_z, source_q, source_w,
                    run_params, pointwisePotential, 0);

        }


    /***************************************/
    /********* Regularized-Yukawa **********/
    /***************************************/

    } else if (run_params->kernel == REGULARIZED_YUKAWA) {

        if (run_params->singularity == SKIPPING) {

            K_RegularizedYukawa_Direct(numTargets, numSources, 0, 0,
                        target_x, target_y, target_z,
                        source_x, source_y, source_z, source_q, source_w,
                        run_params, pointwisePotential, 0);

        } else if (run_params->singularity == SUBTRACTION) {

            K_RegularizedYukawa_SS_Direct(numTargets, numSources, 0, 0,
                        target_x, target_y, target_z, target_q,
                        source_x, source_y, source_z, source_q, source_w,
                        run_params, pointwisePotential, 0);
        }



    /***************************************/
    /********* Atan ************************/
    /***************************************/

    } else if (run_params->kernel == ATAN) {

            K_Atan_Direct(numTargets, numSources, 0, 0,
                        target_x, target_y, target_z,
                        source_x, source_y, source_z, source_q, source_w,
                        run_params, pointwisePotential, 0);

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




void InteractionCompute_SubtractionPotentialCorrection(double *pointwisePotential, double *target_q, int numTargets,
                                  struct RunParams *run_params)
{

    if (run_params->singularity == SUBTRACTION) {
        if (run_params->kernel == COULOMB) {
            K_Coulomb_SS_Correction(pointwisePotential, target_q, numTargets, run_params);

        } else if (run_params->kernel == REGULARIZED_COULOMB) {
            K_RegularizedCoulomb_SS_Correction(pointwisePotential, target_q, numTargets, run_params);

        } else if (run_params->kernel == YUKAWA) {
            K_Yukawa_SS_Correction(pointwisePotential, target_q, numTargets, run_params);

        } else if (run_params->kernel == REGULARIZED_YUKAWA) {
            K_RegularizedYukawa_SS_Correction(pointwisePotential, target_q, numTargets, run_params);

        }
    }

    return;
}
