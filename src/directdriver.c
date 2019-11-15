/*
 *Procedures for Particle-Cluster Treecode
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>

#include "array.h"
#include "globvars.h"
#include "tools.h"

#include "struct_nodes.h"
#include "struct_particles.h"

#include "kernels/coulomb.h"
#include "kernels/yukawa.h"
#include "kernels/coulomb_singularity_subtraction.h"
#include "kernels/yukawa_singularity_subtraction.h"

#include "directdriver.h"


void directdriver(struct particles *sources, struct particles *targets,
                  char *kernelName, double kernel_parameter, char *singularityHandling,
                  char *approximationName, double *pointwisePotential, double *totalPotential)
{
    int rank, numProcs, ierr;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    int numSources = sources->num;
    int numTargets = targets->num;

    double *source_x = sources->x;
    double *source_y = sources->y;
    double *source_z = sources->z;
    double *source_q = sources->q;
    double *source_w = sources->w;

    double *target_x = targets->x;
    double *target_y = targets->y;
    double *target_z = targets->z;
    double *target_q = targets->q;


#ifdef OPENACC_ENABLED
    #pragma acc data copyin(source_x[0:numSources], source_y[0:numSources], source_z[0:numSources], \
                            source_q[0:numSources], source_w[0:numSources], \
                            target_x[0:numTargets], target_y[0:numTargets], target_z[0:numTargets], \
                            target_q[0:numTargets]), copy(pointwisePotential[0:numTargets])
#endif
    {


/**********************************************************/
/************** POTENTIAL FROM DIRECT *********************/
/**********************************************************/


    /***************************************/
    /********* Coulomb *********************/
    /***************************************/

    if (strcmp(kernelName, "coulomb") == 0) {

        if (strcmp(singularityHandling, "skipping") == 0) {

            coulombDirect(numTargets, numSources, 0, 0,
                    target_x, target_y, target_z,
                    source_x, source_y, source_z, source_q, source_w,
                    pointwisePotential, 0);

        } else if (strcmp(singularityHandling, "subtraction") == 0) {

            coulombSingularitySubtractionDirect(numTargets, numSources, 0, 0,
                    target_x, target_y, target_z, target_q,
                    source_x, source_y, source_z, source_q, source_w,
                    kernel_parameter, pointwisePotential, 0);

        }else {
            printf("Invalid choice of singularityHandling. Exiting. \n");
            exit(1);
        }

    /***************************************/
    /********* Yukawa **********************/
    /***************************************/

    } else if (strcmp(kernelName, "yukawa") == 0) {

        if (strcmp(singularityHandling, "skipping") == 0) {

            yukawaDirect(numTargets, numSources, 0, 0,
                    target_x, target_y, target_z,
                    source_x, source_y, source_z, source_q, source_w,
                    kernel_parameter, pointwisePotential, 0);

        } else if (strcmp(singularityHandling, "subtraction") == 0) {

            yukawaSingularitySubtractionDirect(numTargets, numSources, 0, 0,
                    target_x, target_y, target_z, target_q,
                    source_x, source_y, source_z, source_q, source_w,
                    kernel_parameter, pointwisePotential, 0);

        } else {
            printf("Invalid choice of singularityHandling. Exiting. \n");
            exit(1);
        }


    } else {
        printf("Invalid kernelName. Exiting.\n");
        exit(1);
    }

#ifdef OPENACC_ENABLED
        #pragma acc wait
#endif
    } // end acc data region


    *totalPotential = sum(pointwisePotential, numTargets);

    return;

} /* END of function pc_treecode */
