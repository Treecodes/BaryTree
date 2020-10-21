#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "tools.h"
#include "../run_params/run_params.h"
#include "../run_params/struct_run_params.h"

#include "advanced_timings.h"


/*----------------------------------------------------------------------------*/
void Timing_Calculate(double time_tree_glob[3][13], double time_tree[13], double total_time_glob[1], double total_time[1])
{
    
    MPI_Reduce(time_tree, &time_tree_glob[0], 13, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(time_tree, &time_tree_glob[1], 13, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(time_tree, &time_tree_glob[2], 13, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    MPI_Reduce(total_time, total_time_glob, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    return;
}


/*----------------------------------------------------------------------------*/
void Timing_Print(double time_tree_glob[3][13], double total_time_glob[1], struct RunParams *run_params)
{
    int rank, numProcs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    
    if (rank == 0) {
    
    /* Printing direct and treecode time calculations: */
    printf("[BaryTree]\n");
    printf("[BaryTree] ");
    printf("Treecode timing summary (all times in seconds)...\n");
    printf("[BaryTree] ");
    printf("                                       Max                           Avg                          Max/Min\n");
    printf("[BaryTree] ");
    printf("|    Treedriver......................  %9.3e s    (100.00%%)  \n",
                 total_time_glob[0]);
    printf("[BaryTree] ");
    printf("|    |\n");
    printf("[BaryTree] ");
    printf("|    |....Build local tree...........  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_tree_glob[1][0],          time_tree_glob[1][0] / total_time_glob[0]*100.0,
                 time_tree_glob[2][0]/numProcs, time_tree_glob[2][0]/numProcs / total_time_glob[0]*100.0,
                 time_tree_glob[1][0]/time_tree_glob[0][0]);
    printf("[BaryTree] ");
    printf("|    |....Build local batches........  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_tree_glob[1][1],          time_tree_glob[1][1] / total_time_glob[0]*100.0,
                 time_tree_glob[2][1]/numProcs, time_tree_glob[2][1]/numProcs / total_time_glob[0]*100.0,
                 time_tree_glob[1][1]/time_tree_glob[0][1]);
    printf("[BaryTree] ");
    printf("|    |....Build local clusters.......  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_tree_glob[1][2],          time_tree_glob[1][2] / total_time_glob[0]*100.0,
                 time_tree_glob[2][2]/numProcs, time_tree_glob[2][2]/numProcs / total_time_glob[0]*100.0,
                 time_tree_glob[1][2]/time_tree_glob[0][2]);

    if (numProcs > 1) {
    printf("[BaryTree] ");
    printf("|    |....Build LET..................  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_tree_glob[1][3],          time_tree_glob[1][3] / total_time_glob[0]*100.0,
                 time_tree_glob[2][3]/numProcs, time_tree_glob[2][3]/numProcs / total_time_glob[0]*100.0,
                 time_tree_glob[1][3]/time_tree_glob[0][3]);
    }

    printf("[BaryTree] ");
    printf("|    |....Build local lists..........  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_tree_glob[1][4],          time_tree_glob[1][4] / total_time_glob[0]*100.0,
                 time_tree_glob[2][4]/numProcs, time_tree_glob[2][4]/numProcs / total_time_glob[0]*100.0,
                 time_tree_glob[1][4]/time_tree_glob[0][4]);
    printf("[BaryTree] ");
    printf("|    |....Compute local..............  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_tree_glob[1][5],          time_tree_glob[1][5] / total_time_glob[0]*100.0,
                 time_tree_glob[2][5]/numProcs, time_tree_glob[2][5]/numProcs / total_time_glob[0]*100.0,
                 time_tree_glob[1][5]/time_tree_glob[0][5]);

    if (numProcs > 1) {
    printf("[BaryTree] ");
    printf("|    |....Build remote lists.........  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_tree_glob[1][6],          time_tree_glob[1][6] / total_time_glob[0]*100.0,
                 time_tree_glob[2][6]/numProcs, time_tree_glob[2][6]/numProcs / total_time_glob[0]*100.0,
                 time_tree_glob[1][6]/time_tree_glob[0][6]);
    printf("[BaryTree] ");
    printf("|    |....Compute remote.............  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_tree_glob[1][7],          time_tree_glob[1][7] / total_time_glob[0]*100.0,
                 time_tree_glob[2][7]/numProcs, time_tree_glob[2][7]/numProcs / total_time_glob[0]*100.0,
                 time_tree_glob[1][7]/time_tree_glob[0][7]);
    }

    if (run_params->compute_type == CLUSTER_PARTICLE || run_params->compute_type == CLUSTER_CLUSTER) {
    printf("[BaryTree] ");
    printf("|    |....Compute cp2................  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_tree_glob[1][8],          time_tree_glob[1][8] / total_time_glob[0]*100.0,
                 time_tree_glob[2][8]/numProcs, time_tree_glob[2][8]/numProcs / total_time_glob[0]*100.0,
                 time_tree_glob[1][8]/time_tree_glob[0][8]);
    }

    printf("[BaryTree] ");
    printf("|    |....Correct potential..........  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_tree_glob[1][9],          time_tree_glob[1][9] / total_time_glob[0]*100.0,
                 time_tree_glob[2][9]/numProcs, time_tree_glob[2][9]/numProcs / total_time_glob[0]*100.0,
                 time_tree_glob[1][9]/time_tree_glob[0][9]);
    printf("[BaryTree] ");
    printf("|    |....Cleanup....................  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_tree_glob[1][10],          time_tree_glob[1][10] / total_time_glob[0]*100.0,
                 time_tree_glob[2][10]/numProcs, time_tree_glob[2][10]/numProcs / total_time_glob[0]*100.0,
                 time_tree_glob[1][10]/time_tree_glob[0][10]);
    printf("[BaryTree]\n");
    
    if (numProcs > 1) {
    printf("[BaryTree] ");
    printf("((   |....Total setup................  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f ))\n",
                 time_tree_glob[1][11],          time_tree_glob[1][11] / total_time_glob[0]*100.0,
                 time_tree_glob[2][11]/numProcs, time_tree_glob[2][11]/numProcs / total_time_glob[0]*100.0,
                 time_tree_glob[1][11]/time_tree_glob[0][11]);
    printf("[BaryTree] ");
    printf("((   |....Build local clusters.......  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f ))\n",
                 time_tree_glob[1][02],          time_tree_glob[1][02] / total_time_glob[0]*100.0,
                 time_tree_glob[2][02]/numProcs, time_tree_glob[2][02]/numProcs / total_time_glob[0]*100.0,
                 time_tree_glob[1][02]/time_tree_glob[0][02]);
    printf("[BaryTree] ");
    printf("((   |....Total compute..............  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f ))\n",
                 time_tree_glob[1][12],          time_tree_glob[1][12] / total_time_glob[0]*100.0,
                 time_tree_glob[2][12]/numProcs, time_tree_glob[2][12]/numProcs / total_time_glob[0]*100.0,
                 time_tree_glob[1][12]/time_tree_glob[0][12]);
    printf("[BaryTree]\n");
    printf("[BaryTree]\n");
    }
    }

    return;
}


