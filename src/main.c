#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <string.h>

#include "array.h"
#include "treedriver.h"
#include "tools.h"
#include "sort.h"


/* The treedriver routine in Fortran */
int main(int argc, char **argv)
{
    int rank, p;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    /* runtime parameters */
    int numparsS, numparsT, order;
    int maxparnode, treelevel;
    int iflag, pot_type, tree_type;

    double theta, temp;
    double kappa;

    /* arrays for coordinates, charges, energy of target particles */
    double *xS, *yS, *zS, *qS;  /* source particles */
    double *xT, *yT, *zT;       /* target particles */
    double *denergy, *tenergy;  /* exact energy, treecode energy */

    /* for potential energy calculation */
    double tpeng, dpeng;
    double dpengglob, tpengglob;

    /* insert variables for date-time calculation? */
    double time_direct, time_tree;
    double time_tree_min, time_tree_max, time_tree_tot;
    double time1, time2;

    /* input and output files */
    char *sampin1, *sampin2, *sampin3, *sampout;

    /* variables for error calculations */
    double inferr, relinferr, n2err, reln2err;
    double inferrglob, relinferrglob, n2errglob, reln2errglob;

    /* local variables */
    int i, j;
    int numparsTloc, maxparsTloc, globparsTloc;
    double buf[5];
    
    /* MPI Variables */
    MPI_File fpmpi;
    MPI_Status status;


    /* Executable statements begin here */

    sampin1 = argv[1];
    if (strcmp(sampin1,"--help") == 0)
    {
        if (rank == 0)
        {
            printf("Input arguments: \n");
            printf("       sampin1:  sources input file \n");               // "S10000.txt"
            printf("       sampin2:  targets input file \n");               // "T1000000.txt"
            printf("       sampin3:  direct calc potential input file \n"); // "ex_s4_t6.txt"
            printf("       sampout:  tree calc potential output file \n");  // "out.txt"
            printf("      numparsS:  number of sources \n");                // 10000
            printf("      numparsT:  number of targets \n");                // 1000000
            printf("         theta:  multipole acceptance criterion \n");   // 0.75
            printf("         order:  order of treecode Taylor expansion \n");        // 20
            printf("     tree_type:  0--cluster-particle, 1--particle-cluster \n");  // 0
            printf("    maxparnode:  maximum particles in leaf \n");                 // 500
            printf("     treelevel:  maximum tree levels \n");                       // 5
            printf("         iflag:  0--use maxparnode, 1--use treelevel \n");       // 0
            printf("         kappa:  screened Coulomb parameter \n");                // 0.00
            printf("      pot_type:  0--Coulomb, 1--screened Coulomb \n");           // 1
        }
        return 0;
    }
    
    sampin2 = argv[2];
    sampin3 = argv[3];
    sampout = argv[4];
    numparsS = atoi(argv[5]);
    numparsT = atoi(argv[6]);
    theta = atof(argv[7]);
    order = atoi(argv[8]);
    tree_type = atoi(argv[9]);
    maxparnode = atoi(argv[10]);
    treelevel = atoi(argv[11]);
    iflag = atoi(argv[12]);
    kappa = atof(argv[13]);
    pot_type = atoi(argv[14]);
    
    
    numparsTloc = (int)floor((double)numparsT/(double)p);
    maxparsTloc = numparsTloc + (numparsT - (int)floor((double)numparsT/(double)p) * p);
    
    if (rank == 0)
        numparsTloc = maxparsTloc;
    
    globparsTloc = maxparsTloc + numparsTloc * (rank-1);


    make_vector(xS,numparsS);
    make_vector(yS,numparsS);
    make_vector(zS,numparsS);
    make_vector(qS,numparsS);
    
    make_vector(xT,numparsTloc);
    make_vector(yT,numparsTloc);
    make_vector(zT,numparsTloc);

    make_vector(tenergy,numparsTloc);
    make_vector(denergy,numparsTloc);


    /* Reading in coordinates and charges for the source particles*/
    MPI_File_open(MPI_COMM_WORLD, sampin1, MPI_MODE_RDONLY, MPI_INFO_NULL, &fpmpi);
    MPI_File_seek(fpmpi, (MPI_Offset)0, MPI_SEEK_SET);
    for (i = 0; i < numparsS; i++) {
        MPI_File_read(fpmpi, buf, 4, MPI_DOUBLE, &status);
        xS[i] = buf[0];
        yS[i] = buf[1];
        zS[i] = buf[2];
        qS[i] = buf[3];
    }
    MPI_File_close(&fpmpi);
    

    /* Reading in coordinates for the targets */
    MPI_File_open(MPI_COMM_WORLD, sampin2, MPI_MODE_RDONLY, MPI_INFO_NULL, &fpmpi);
    MPI_File_seek(fpmpi, (MPI_Offset)globparsTloc*3*sizeof(double), MPI_SEEK_SET);
    for (i = 0; i < numparsTloc; i++) {
        MPI_File_read(fpmpi, buf, 3, MPI_DOUBLE, &status);
        xT[i] = buf[0];
        yT[i] = buf[1];
        zT[i] = buf[2];
    }
    MPI_File_close(&fpmpi);

    
    /* Reading in the values for exact energy at each target */
    MPI_File_open(MPI_COMM_WORLD, sampin3, MPI_MODE_RDONLY, MPI_INFO_NULL, &fpmpi);
    if (rank == 0)
    {
        MPI_File_seek(fpmpi, (MPI_Offset)0, MPI_SEEK_SET);
        MPI_File_read(fpmpi, &time_direct, 1, MPI_DOUBLE, &status);
    }
    MPI_File_seek(fpmpi, (MPI_Offset)(globparsTloc+1)*sizeof(double), MPI_SEEK_SET);
    for (i = 0; i < numparsTloc; i++) {
        MPI_File_read(fpmpi, buf, 1, MPI_DOUBLE, &status);
        denergy[i] = buf[0];
    }
    MPI_File_close(&fpmpi);
    
    time1 = MPI_Wtime();
    sortTargets(xT, yT, zT, numparsTloc);
    time2 = MPI_Wtime();
    printf("  Sort time (s): %f\n\n", time2-time1);
  
    
    
    /* Calling main treecode subroutine to calculate approximate energy */
    treecode(xS, yS, zS, qS, xT, yT, zT, numparsS, numparsTloc,
             tenergy, &tpeng, order, theta, 1, maxparnode, &time_tree,
             treelevel, iflag, pot_type, kappa, tree_type);

    
    
    /* Reducing values to root process */
    MPI_Reduce(&time_tree, &time_tree_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&time_tree, &time_tree_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&time_tree, &time_tree_tot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    dpeng = sum(denergy, numparsTloc);
    MPI_Reduce(&dpeng, &dpengglob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&tpeng, &tpengglob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    
    if (rank == 0)
    {
        /* Printing direct and treecode time calculations: */
        printf("                   Direct time (s):  %f\n", time_direct);
        printf("                 Max tree time (s):  %f\n", time_tree_max);
        printf("                 Min tree time (s):  %f\n", time_tree_min);
        printf("                 Avg tree time (s):  %f\n\n", time_tree_tot/(double)p);
        printf("         Direct : Tree on %d procs:  %f\n\n",
               p, time_direct/(time_tree_max*(double)p));

    
        /* Printing error in potential energy and potential energy */
        printf("           Direct potential energy:  %f\n", dpengglob);
        printf("             Tree potential energy:  %f\n\n", tpengglob);
    
        printf("Absolute error for total potential:  %e\n",
               fabs(tpengglob-dpengglob));
        printf("Relative error for total potential:  %e\n\n",
               fabs((tpengglob-dpengglob)/dpengglob));
    }
    
    /* Computing pointwise potential errors */

    inferr = 0.0;
    relinferr = 0.0;
    n2err = 0.0;
    reln2err = 0.0;

    for (j = 0; j < numparsTloc; j++) {
        temp = fabs(denergy[j] - tenergy[j]);
        
        if (temp >= inferr)
            inferr = temp;

        if (fabs(denergy[j]) >= relinferr)
            relinferr = fabs(denergy[j]);

        n2err = n2err + pow(denergy[j] - tenergy[j], 2.0);
        reln2err = reln2err + pow(denergy[j], 2.0);
    }

    MPI_Reduce(&inferr, &inferrglob, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&relinferr, &relinferrglob, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    
    MPI_Reduce(&n2err, &n2errglob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&reln2err, &reln2errglob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    
    if (rank == 0)
    {
        relinferrglob = inferrglob / relinferrglob;
        reln2errglob = sqrt(n2errglob / reln2errglob);
        n2errglob = sqrt(n2errglob);

        printf("Absolute inf norm error in potential:  %e \n", inferrglob);
        printf("Relative inf norm error in potential:  %e \n\n", relinferrglob);
        printf("  Absolute 2 norm error in potential:  %e \n", n2errglob);
        printf("  Relative 2 norm error in potential:  %e \n\n", reln2errglob);
    }
    
    MPI_Finalize();
    return 0;
}
