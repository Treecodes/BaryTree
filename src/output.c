
/*
 *Procedures for Output Struct
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>

#include "array.h"
#include "globvars.h"
#include "struct_nodes.h"
#include "struct_particles.h"
#include "struct_clusters.h"
#include "tools.h"

#include "output.h"



void Output_Alloc(struct output *output, int numberOfParticles, int forces)
{
    output->numberOfParticles = numberOfParticles;
    make_vector(output->potential, output->numberOfParticles);

    if (forces==1){
        make_vector(output->forcesX, output->numberOfParticles);
        make_vector(output->forcesY, output->numberOfParticles);
        make_vector(output->forcesZ, output->numberOfParticles);
    }

    return;
}   /* END of function Output_Alloc */




void Output_Free(struct output *output)
{
    free_vector(output->potential);
    if (output->forcesX != NULL) {
        free_vector(output->forcesX);
        free_vector(output->forcesY);
        free_vector(output->forcesZ);
    }
    free(output);

    return;
}   /* END of function Output_Free */
