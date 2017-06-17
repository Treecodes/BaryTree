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
    int maxparnode, pot_type;
    int sflag, dflag;

    double theta, temp;
    double kappa;

    /* source particles */
    double *xS = NULL;
    double *yS = NULL;
    double *zS = NULL;
    double *qS = NULL;
    
    /* target particles */
    double xyzminmax[6];
    int xyzdim[3];
    
    /* exact energy, treecode energy */
    double *denergy = NULL;
    double *tenergy = NULL;

    /* for potential energy calculation */
    double tpeng = 0;
    double dpeng = 0;
    double tpengglob = 0;
    double dpengglob = 0;

    /* insert variables for date-time calculation? */
    double time_direct, time_tree[4], time_preproc;
    double time_tree_glob[3][4];
    double time1, time2;

    /* input and output files */
    char *sampin1 = NULL;
    char *sampin3 = NULL;
    char *sampout = NULL;
    FILE *fp;

    /* variables for error calculations */
    double inferr, relinferr, n2err, reln2err;

    /* local variables */
    int i, j;
    double buf[5];
    
    int numparsTloc;
    int xyzdimloc[3];
    double xyzminmaxloc[3];
    
    int *displs = NULL;
    int *scounts = NULL;
    
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
            printf("       sampin3:  direct calc potential input file \n"); // "ex_s4_t6.txt"
            printf("       sampout:  tree calc potential output file \n");  // "out.txt"
            printf("      numparsS:  number of sources \n");                // 10000
            printf("         theta:  multipole acceptance criterion \n");   // 0.75
            printf("         order:  order of treecode Taylor expansion \n");    // 20
            printf("    maxparnode:  maximum particles in leaf \n");             // 500
            printf("         kappa:  screened Coulomb parameter \n");            // 0.00
            printf("      pot_type:  0--Coulomb, 1--screened Coulomb \n");       // 1

            printf("         sflag:  0--sort, 1--interleave \n");                // 0
            printf("         dflag:  direction 0--x, 1--y, 2--z \n");            // 0
            
            printf("          xmin:  if on grid, min x dimension \n");           // 0
            printf("          xmax:  if on grid, max x dimension \n");           // 0
            printf("          ymin:  if on grid, min y dimension \n");           // 0
            printf("          ymax:  if on grid, max y dimension \n");           // 0
            printf("          zmin:  if on grid, min z dimension \n");           // 0
            printf("          zmax:  if on grid, max z dimension \n");           // 0
            
            printf("          xdim:  if on grid, number x gridpoints \n");       // 0
            printf("          ydim:  if on grid, number y gridpoints \n");       // 0
            printf("          zdim:  if on grid, number z gridpoints \n");       // 0
        }
        return 0;
    }
    
    sampin3 = argv[2];
    sampout = argv[3];
    numparsS = atoi(argv[4]);
    theta = atof(argv[5]);
    order = atoi(argv[6]);
    maxparnode = atoi(argv[7]);
    kappa = atof(argv[8]);
    pot_type = atoi(argv[9]);
    
    sflag = atoi(argv[10]);
    dflag = atoi(argv[11]);

    xyzminmax[0] = atof(argv[12]);
    xyzminmax[1] = atof(argv[13]);
    xyzminmax[2] = atof(argv[14]);
    xyzminmax[3] = atof(argv[15]);
    xyzminmax[4] = atof(argv[16]);
    xyzminmax[5] = atof(argv[17]);

    xyzdim[0] = atoi(argv[18]);
    xyzdim[1] = atoi(argv[19]);
    xyzdim[2] = atoi(argv[20]);

    
    time1 = MPI_Wtime();
    numparsT = xyzdim[0] * xyzdim[1] * xyzdim[2];
    
    
    
    
    
    make_vector(xS,numparsS);
    make_vector(yS,numparsS);
    make_vector(zS,numparsS);
    make_vector(qS,numparsS);
    
    make_vector(tenergy,numparsT);
    make_vector(denergy,numparsT);
        
    
    MPI_File_open(MPI_COMM_SELF, sampin3, MPI_MODE_RDONLY, MPI_INFO_NULL, &fpmpi);
    MPI_File_seek(fpmpi, (MPI_Offset)0, MPI_SEEK_SET);
    MPI_File_read(fpmpi, &time_direct, 1, MPI_DOUBLE, &status);
    MPI_File_read(fpmpi, denergy, numparsT, MPI_DOUBLE, &status);
    MPI_File_close(&fpmpi);
    

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

    time2 = MPI_Wtime();
    time_preproc = time2 - time1;


    /* Calling main treecode subroutine to calculate approximate energy */
    treecode_grid(xS, yS, zS, qS, xyzminmax, xyzdim, numparsS,
                  tenergy, &tpengglob, order, theta, maxparnode, time_tree,
                  pot_type, kappa);

    dpengglob = sum(denergy, numparsT);
    
    
    if (rank == 0) {
        /* Printing direct and treecode time calculations: */
        printf("                   Direct time (s):  %f\n\n", time_direct);
        printf("              Pre-process time (s):  %f\n", time_preproc);
        printf("      Min, Max tree setup time (s):  %f\n", time_tree[0]);
        printf("             Min, Max cp1 time (s):  %f\n", time_tree[1]);
        printf("             Min, Max cp2 time (s):  %f\n", time_tree[2]);
        
        printf("      Min, Max total tree time (s):  %f\n\n", time_tree[3]);
        printf(" Preproc + Max total tree time (s):  %f\n\n", time_tree[3] + time_preproc);
    
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

    for (j = 0; j < numparsT; j++) {
        temp = fabs(denergy[j] - tenergy[j]);
    
        if (temp >= inferr) inferr = temp;
        if (fabs(denergy[j]) >= relinferr) relinferr = fabs(denergy[j]);

        n2err = n2err + pow(denergy[j] - tenergy[j], 2.0);
        reln2err = reln2err + pow(denergy[j], 2.0);
    }

    relinferr = inferr / relinferr;
    reln2err = sqrt(n2err / reln2err);
    n2err = sqrt(n2err);

    
    if (rank == 0) {
        printf("Absolute inf norm error in potential:  %e \n", inferr);
        printf("Relative inf norm error in potential:  %e \n\n", relinferr);
        printf("  Absolute 2 norm error in potential:  %e \n", n2err);
        printf("  Relative 2 norm error in potential:  %e \n\n", reln2err);
    
        fp = fopen(sampout, "a");
        fprintf(fp, "%s \t %s \t %d \t %f \t %f \t %f \t %f \t %f \t"
                    "%f \t %d \t %d \t %d \t %f \t %d \t"
                    "%d \t %f \t %d \t %f \t %f \t %f \t %f \t %f \t"
                    "%f \t %f \t %f \t %f \t %e \t %e \t %e \t %e \n",
                sampin1, sampin3, numparsS, xyzminmax[0], xyzminmax[1], xyzminmax[2],
                xyzminmax[3], xyzminmax[4], xyzminmax[5], xyzdim[0], xyzdim[1], xyzdim[2],
                theta, order, maxparnode, kappa, pot_type,
                time_preproc, time_tree[0], time_tree[1], time_tree[2], time_tree[3],
                dpengglob, tpengglob, fabs(tpengglob-dpengglob),
                fabs((tpengglob-dpengglob)/dpengglob),
                inferr, relinferr, n2err, reln2err);
        fclose(fp);
    }
    
    
    free_vector(xS);
    free_vector(yS);
    free_vector(zS);
    free_vector(qS);
    
    free_vector(denergy);
    free_vector(tenergy);
    
    MPI_Finalize();
    return 0;
    
}
