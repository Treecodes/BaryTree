#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <string.h>

#include "array.h"
#include "tools.h"
#include "treedriver.h"

#include "struct_particles.h"


int main(int argc, char **argv)
{
    int rank, numProcs;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    /* runtime parameters */
    int verbosity=0;
    int numparsS, numparsT, order;
    int maxparnode, batch_size;
    double theta, kappa;

    struct particles *sources = NULL;
    struct particles *targets = NULL;

    double *tenergy = NULL;
    double *denergy = NULL;

    /* for potential energy calculation */
    double dpeng = 0, tpeng = 0;
    double dpengglob = 0, tpengglob = 0;

    /* variables for date-time calculation */
    double time_run[3], time_tree[10], time_direct;
    double time_run_glob[3][3], time_tree_glob[3][10];
    double time1, time2;

    /* input and output files */
    char *sampin1 = NULL;
    char *sampin2 = NULL;
    char *sampin3 = NULL;
    char *offset1 = NULL;
    char *offset2 = NULL;
    char *sampout = NULL;
    char *kernelName = NULL;
    char *singularityHandling = NULL;
    char *approximationName = NULL;
    FILE *fp;

    double buf[5];
    
    double *xS = NULL;
    double *yS = NULL;
    double *zS = NULL;
    double *qS = NULL;
    double *wS = NULL;

    double *xT = NULL;
    double *yT = NULL;
    double *zT = NULL;
    double *qT = NULL;

    /* MPI Variables */
    int *displs = NULL;
    int *scounts = NULL;
    int mpi_err;
    MPI_File fpmpi;
    MPI_Status status;


    /* Executable statements begin here */
    sampin1 = argv[1];
    if (strcmp(sampin1,"--help") == 0) {
        if (rank == 0) {
            printf("Input arguments... \n");
            printf("        infile 1:  sources input file \n");
            printf("        infile 2:  targets input file \n");
            printf("        infile 3:  sources offset file \n");
            printf("        infile 4:  targets offset file \n");
            printf("        infile 5:  direct calc potential input file \n");
            printf("     csv outfile:  results summary to CSV file \n");
            printf("        numparsS:  number of sources \n");
            printf("        numparsT:  number of targets \n");
            printf("           theta:  multipole acceptance criterion \n");
            printf("           order:  number of Chebyshev interp. pts per Cartesian direction \n");
            printf("      maxparnode:  maximum particles in leaf \n");
            printf("      batch size:  maximum size of target batch \n");
            printf(" pot/approx type:  0--Coulomb, Lagrange approx.\n");
            printf("                   1--screened Coulomb/Yukawa, Lagrange approx.\n");
            printf("                   4--Coulomb, Hermite approx.\n");
            printf("                   5--screened Coulomb/Yukawa, Hermite approx.\n");
            printf("           kappa:  screened Coulomb parameter \n");
        }
        return 0;
    }
    
    sampin2 = argv[2];
    offset1 = argv[3];
    offset2 = argv[4];
    sampin3 = argv[5];
    sampout = argv[6];
    numparsS = atoi(argv[7]);
    numparsT = atoi(argv[8]);
    theta = atof(argv[9]);
    order = atoi(argv[10]);
    maxparnode = atoi(argv[11]);
    batch_size = atoi(argv[12]);
    kappa = atof(argv[14]);
    
    kernelName = argv[13];
    singularityHandling = argv[15];

    printf("singularityHandling = %s\n", singularityHandling);





    int numparsSloc, numparsTloc, local_sources_offset, local_targets_offset;
    
    if ((mpi_err = MPI_File_open(MPI_COMM_WORLD, offset1,
        MPI_MODE_RDONLY, MPI_INFO_NULL, &fpmpi)) != MPI_SUCCESS) {
        numparsSloc = (int)floor((double)numparsS/(double)numProcs);
        local_sources_offset = numparsSloc * rank;
        if (rank == numProcs-1)
            numparsSloc += (numparsS - (int)floor((double)numparsS/(double)numProcs) * numProcs);
    } else {
        int *sources_offset;
        make_vector(sources_offset, numProcs);
        MPI_File_read(fpmpi, sources_offset, numProcs, MPI_INT, &status);
        if (rank < numProcs-1)
            numparsSloc = sources_offset[rank+1] - sources_offset[rank];
        else
            numparsSloc = numparsS - sources_offset[rank];
        local_sources_offset = sources_offset[rank];
        MPI_File_close(&fpmpi);
        free_vector(sources_offset);
    }
    
    if ((mpi_err = MPI_File_open(MPI_COMM_WORLD, offset1,
        MPI_MODE_RDONLY, MPI_INFO_NULL, &fpmpi)) != MPI_SUCCESS) {
        numparsTloc = (int)floor((double)numparsT/(double)numProcs);
        local_targets_offset = numparsTloc * rank;
        if (rank == numProcs-1)
            numparsTloc += (numparsT - (int)floor((double)numparsT/(double)numProcs) * numProcs);
    } else {
        int *targets_offset;
        make_vector(targets_offset, numProcs);
        MPI_File_read(fpmpi, targets_offset, numProcs, MPI_INT, &status);
        if (rank < numProcs-1)
            numparsTloc = targets_offset[rank+1] - targets_offset[rank];
        else
            numparsTloc = numparsT - targets_offset[rank];
        local_targets_offset = targets_offset[rank];
        MPI_File_close(&fpmpi);
        free_vector(targets_offset);
    }

    
    time1 = MPI_Wtime();
    
    sources = malloc(sizeof(struct particles));
    targets = malloc(sizeof(struct particles));

    make_vector(tenergy, numparsTloc);
    make_vector(denergy, numparsTloc);

    sources->num = numparsSloc;
    make_vector(sources->x, numparsSloc);
    make_vector(sources->y, numparsSloc);
    make_vector(sources->z, numparsSloc);
    make_vector(sources->q, numparsSloc);
    make_vector(sources->w, numparsSloc);

    targets->num = numparsTloc;
    make_vector(targets->x, numparsTloc);
    make_vector(targets->y, numparsTloc);
    make_vector(targets->z, numparsTloc);
    make_vector(targets->q, numparsTloc);
    make_vector(targets->order, numparsTloc);
    make_vector(tenergy, numparsTloc);

    /* Reading in coordinates and charges for the source particles*/
    if ((mpi_err = MPI_File_open(MPI_COMM_WORLD, sampin1, MPI_MODE_RDONLY,
                  MPI_INFO_NULL, &fpmpi)) != MPI_SUCCESS) {
        printf("Error! Could not open sources input file. Exiting.\n");
        return 1;
    }
    
    MPI_File_seek(fpmpi, (MPI_Offset) (local_sources_offset * 5 * sizeof(double)), MPI_SEEK_SET);
    for (int i = 0; i < numparsSloc; ++i) {
        MPI_File_read(fpmpi, buf, 5, MPI_DOUBLE, &status);
        sources->x[i] = buf[0];
        sources->y[i] = buf[1];
    sources->z[i] = buf[2];
    sources->q[i] = buf[3];
    sources->w[i] = buf[4];
    }
    MPI_File_close(&fpmpi);

    /* Reading in coordinates for target particles*/
    if ((mpi_err = MPI_File_open(MPI_COMM_WORLD, sampin2, MPI_MODE_RDONLY,
                  MPI_INFO_NULL, &fpmpi)) != MPI_SUCCESS) {
        printf("Error! Could not open targets input file. Exiting.\n");
        return 1;
    }
    
    MPI_File_seek(fpmpi, (MPI_Offset) (local_targets_offset * 4 * sizeof(double)), MPI_SEEK_SET);
    for (int i = 0; i < numparsTloc; ++i) {
        MPI_File_read(fpmpi, buf, 4, MPI_DOUBLE, &status);
        targets->x[i] = buf[0];
        targets->y[i] = buf[1];
        targets->z[i] = buf[2];
        targets->q[i] = buf[3];
        targets->order[i] = i;
    }
    MPI_File_close(&fpmpi);



    // reading in file containing direct sum results.
    MPI_File_open(MPI_COMM_SELF, sampin3, MPI_MODE_RDONLY, MPI_INFO_NULL, &fpmpi);
    if (rank == 0) MPI_File_read(fpmpi, &time_direct, 1, MPI_DOUBLE, &status);
    MPI_File_seek(fpmpi, (MPI_Offset)((1 + local_targets_offset) * sizeof(double)), MPI_SEEK_SET);
    MPI_File_read(fpmpi, denergy, numparsTloc, MPI_DOUBLE, &status);
    MPI_File_close(&fpmpi);


    // Set up kernel
    for (int i = 0; i < numparsTloc; i++) tenergy[i] = 0.0;


#ifdef OPENACC_ENABLED
    #pragma acc set device_num(rank) device_type(acc_device_nvidia)
    #pragma acc init device_type(acc_device_nvidia)
#endif

    time_run[0] = MPI_Wtime() - time1;

    /* Calling main treecode subroutine to calculate approximate energy */

    MPI_Barrier(MPI_COMM_WORLD);
    
    time1 = MPI_Wtime();
    
    treedriver(sources, targets, order, theta, maxparnode, batch_size,
               kernelName, kappa, singularityHandling, approximationName, 1, tenergy, time_tree, verbosity);

    tpeng = sum(tenergy, targets->num);
               
    time_run[1] = MPI_Wtime() - time1;
    time_run[2] = time_run[0] + time_run[1];
    
    MPI_Barrier(MPI_COMM_WORLD);

    
    /* Reducing values to root process */
    MPI_Reduce(time_tree, &time_tree_glob[0], 10, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(time_tree, &time_tree_glob[1], 10, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(time_tree, &time_tree_glob[2], 10, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    MPI_Reduce(time_run, &time_run_glob[0], 3, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(time_run, &time_run_glob[1], 3, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(time_run, &time_run_glob[2], 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    dpeng = sum(denergy, numparsTloc);
    MPI_Barrier(MPI_COMM_WORLD);
    
    MPI_Reduce(&tpeng, &tpengglob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&dpeng, &dpengglob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    
    if (rank == 0)
    {
    
        double min_percent = 100. / time_run_glob[0][2];
        double max_percent = 100. / time_run_glob[1][2];
        double avg_percent = 100. / time_run_glob[2][2];

        /* Printing direct and treecode time calculations: */
        printf("\n\nTreecode timing summary (all times in seconds)...\n\n");
        printf("                                       Avg                           Min                         Max\n");
        printf("|    Total time......................  %9.3e s    (100.00%%)      %9.3e s    (100.00%%)    %9.3e s    (100.00%%) \n",
                     time_run_glob[2][2]/numProcs, time_run_glob[0][2], time_run_glob[1][2]);
        printf("|    |\n");
        printf("|    |....Pre-process................  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %9.3e s    (%6.2f%%) \n",
                     time_run_glob[2][0]/numProcs, time_run_glob[2][0] * avg_percent,
                     time_run_glob[0][0], time_run_glob[0][0] * min_percent,
                     time_run_glob[1][0], time_run_glob[1][0] * max_percent);
        printf("|    |....Treedriver.................  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %9.3e s    (%6.2f%%) \n",
                     time_run_glob[2][1]/numProcs, time_run_glob[2][1] * avg_percent,
                     time_run_glob[0][1], time_run_glob[0][1] * min_percent,
                     time_run_glob[1][1], time_run_glob[1][1] * max_percent);
        printf("|         |\n");
        printf("|         |....Build local tree......  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %9.3e s    (%6.2f%%) \n",
                     time_tree_glob[2][0]/numProcs, time_tree_glob[2][0] * avg_percent,
                     time_tree_glob[0][0],          time_tree_glob[0][0] * min_percent,
                     time_tree_glob[1][0],          time_tree_glob[1][0] * max_percent);
        printf("|         |....Build local batches...  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %9.3e s    (%6.2f%%) \n",
                     time_tree_glob[2][1]/numProcs, time_tree_glob[2][1] * avg_percent,
                     time_tree_glob[0][1],          time_tree_glob[0][1] * min_percent,
                     time_tree_glob[1][1],          time_tree_glob[1][1] * max_percent);
        printf("|         |....Fill local clusters...  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %9.3e s    (%6.2f%%) \n",
                     time_tree_glob[2][2]/numProcs, time_tree_glob[2][2] * avg_percent,
                     time_tree_glob[0][2],          time_tree_glob[0][2] * min_percent,
                     time_tree_glob[1][2],          time_tree_glob[1][2] * max_percent);

        if (numProcs > 1) {
        printf("|         |....Build LET.............  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %9.3e s    (%6.2f%%) \n",
                     time_tree_glob[2][3]/numProcs, time_tree_glob[2][3] * avg_percent,
                     time_tree_glob[0][3],          time_tree_glob[0][3] * min_percent,
                     time_tree_glob[1][3],          time_tree_glob[1][3] * max_percent);
        }

        printf("|         |....Build local lists.....  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %9.3e s    (%6.2f%%) \n",
                     time_tree_glob[2][4]/numProcs, time_tree_glob[2][4] * avg_percent,
                     time_tree_glob[0][4],          time_tree_glob[0][4] * min_percent,
                     time_tree_glob[1][4],          time_tree_glob[1][4] * max_percent);
        printf("|         |....Compute local.........  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %9.3e s    (%6.2f%%) \n",
                     time_tree_glob[2][5]/numProcs, time_tree_glob[2][5] * avg_percent,
                     time_tree_glob[0][5],          time_tree_glob[0][5] * min_percent,
                     time_tree_glob[1][5],          time_tree_glob[1][5] * max_percent);

        if (numProcs > 1) {
        printf("|         |....Build remote lists....  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %9.3e s    (%6.2f%%) \n",
                     time_tree_glob[2][6]/numProcs, time_tree_glob[2][6] * avg_percent,
                     time_tree_glob[0][6],          time_tree_glob[0][6] * min_percent,
                     time_tree_glob[1][6],          time_tree_glob[1][6] * max_percent);
        printf("|         |....Compute remote........  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %9.3e s    (%6.2f%%) \n",
                     time_tree_glob[2][7]/numProcs, time_tree_glob[2][7] * avg_percent,
                     time_tree_glob[0][7],          time_tree_glob[0][7] * min_percent,
                     time_tree_glob[1][7],          time_tree_glob[1][7] * max_percent);
        }

        printf("|         |....Correct potential.....  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %9.3e s    (%6.2f%%) \n",
                     time_tree_glob[2][8]/numProcs, time_tree_glob[2][8] * avg_percent,
                     time_tree_glob[0][8],          time_tree_glob[0][8] * min_percent,
                     time_tree_glob[1][8],          time_tree_glob[1][8] * max_percent);
        printf("|         |....Cleanup...............  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %9.3e s    (%6.2f%%) \n\n",
                     time_tree_glob[2][9]/numProcs, time_tree_glob[2][9] * avg_percent,
                     time_tree_glob[0][9],          time_tree_glob[0][9] * min_percent,
                     time_tree_glob[1][9],          time_tree_glob[1][9] * max_percent);
        
        printf("\nTreecode comparison summary...\n\n");
        printf("   Speedup over reported direct time:  %fx\n\n", time_direct/time_run_glob[1][2]);
        printf("             Direct potential energy:  %f\n", dpengglob);
        printf("               Tree potential energy:  %f\n\n", tpengglob);
    
        printf("  Absolute error for total potential:  %e\n",
               fabs(tpengglob-dpengglob));
        printf("  Relative error for total potential:  %e\n\n",
               fabs((tpengglob-dpengglob)/dpengglob));
    }
    
    
    double glob_reln2_err, glob_relinf_err, glob_n2_err, glob_inf_err;
    double inferr = 0.0, relinferr = 0.0, n2err = 0.0, reln2err = 0.0;
    double x, y, z, temp;

    for (int j = 0; j < numparsTloc; ++j) {
        temp = fabs(denergy[targets->order[j]] - tenergy[j]);

        if (temp >= inferr) {
            inferr = temp;
            x = targets->x[targets->order[j]];
            y = targets->y[targets->order[j]];
            z = targets->z[targets->order[j]];
        }

        if (fabs(denergy[j]) >= relinferr)
            relinferr = fabs(denergy[j]);

        n2err = n2err + pow(denergy[targets->order[j]]
                          - tenergy[j], 2.0)*sources->w[j];
        reln2err = reln2err + pow(denergy[j], 2.0)*sources->w[j];
    }

    MPI_Reduce(&reln2err, &glob_reln2_err, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&n2err, &glob_n2_err, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&relinferr, &glob_relinf_err, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&inferr, &glob_inf_err, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        glob_reln2_err = sqrt(glob_n2_err / glob_reln2_err);
        glob_n2_err = sqrt(glob_n2_err);
        glob_relinf_err = glob_inf_err / glob_relinf_err;
        printf("Relative inf norm error in potential:  %e \n", glob_relinf_err);
        printf("  Relative 2 norm error in potential:  %e \n\n", glob_reln2_err);
    }
    
    if (rank == 0) {
        fp = fopen(sampout, "a");
        fprintf(fp, "%s,%s,%s,%d,%d,%f,%d,%d,%d,%f,%s,%d,"
                    "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,"
                    "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,"
                    "%e,%e,%e,%e,%e,%e,%e,%e\n",
            sampin1, sampin2, sampin3, numparsS, numparsT, theta, order,
            maxparnode, batch_size, kappa, kernelName, numProcs, // 1 ends
            time_run_glob[0][0],  time_run_glob[1][0],  // min, max, avg pre-process
            time_run_glob[2][0]/numProcs,
            time_run_glob[0][1],  time_run_glob[1][1],  // min, max, avg treedriver
            time_run_glob[2][1]/numProcs,
            time_tree_glob[0][0], time_tree_glob[1][0],     // min, max, avg build local tree
            time_tree_glob[2][0]/numProcs,
            time_tree_glob[0][1], time_tree_glob[1][1],     // min, max, avg build local batches
            time_tree_glob[2][1]/numProcs,
            time_tree_glob[0][2], time_tree_glob[1][2],     // min, max, avg fill local clusters
            time_tree_glob[2][2]/numProcs,
            time_tree_glob[0][3], time_tree_glob[1][3],     // min, max, avg build LET
            time_tree_glob[2][3]/numProcs,
            time_tree_glob[0][4], time_tree_glob[1][4],     // min, max, avg build local lists
            time_tree_glob[2][4]/numProcs, // 2 ends
            time_tree_glob[0][5], time_tree_glob[1][5],     // min, max, avg compute local interations 
            time_tree_glob[2][5]/numProcs,
            time_tree_glob[0][6], time_tree_glob[1][6],     // min, max, avg build remote lists
            time_tree_glob[2][6]/numProcs,
            time_tree_glob[0][7], time_tree_glob[1][7],     // min, max, avg compute remote interactions
            time_tree_glob[2][7]/numProcs,
            time_tree_glob[0][8], time_tree_glob[1][8],     // min, max, avg correct potential
            time_tree_glob[2][8]/numProcs,
            time_tree_glob[0][9], time_tree_glob[1][9],     // min, max, avg cleanup
            time_tree_glob[2][9]/numProcs,
            time_run_glob[0][2],  time_run_glob[1][2],  // min, max, avg total time
            time_run_glob[2][2]/numProcs, // 3 ends
            dpengglob, tpengglob, fabs(tpengglob-dpengglob),
            fabs((tpengglob-dpengglob)/dpengglob),
            glob_inf_err, glob_relinf_err, glob_n2_err, glob_reln2_err); // 4 ends
        fclose(fp);
    }
    
    free_vector(sources->x);
    free_vector(sources->y);
    free_vector(sources->z);
    free_vector(sources->q);
    free_vector(sources->w);

    free_vector(targets->x);
    free_vector(targets->y);
    free_vector(targets->z);
    free_vector(targets->q);
    free_vector(targets->order);
    
    free(sources);
    free(targets);
    free_vector(tenergy);

    MPI_Finalize();

    return 0;
}
