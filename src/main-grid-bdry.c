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
    double xyzminmax[6], xxminmax[3][4];
    double xyzdd[3];
    int xyzdim[3], xxdim[3][2];
    
    /* exact energy, treecode energy */
    double *denergyglob = NULL;
    double *tenergyglob = NULL;
    double *tenergy = NULL;

    /* for potential energy calculation */
    double tpeng = 0;
    double tpengg[6];
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
    int i, j, tempdim;
    double buf[5];
    
    int numparsTloc, numparsSloc;
    int xyzdimloc[3];
    double xyzddloc[3];
    double xyzminmaxloc[6];
    
    int *displs = NULL;
    int *scounts = NULL;
    int *displsl = NULL;
    int *scountsl = NULL;
    
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

            printf("         sflag:  0--target sort, 1--target interleave, 2--sources \n");  // 0
            printf("                 (note: sflag for determining parallel distribution.\n");
            printf("                        currently only sources implemented.)\n");
            printf("         dflag:  if targets, direction 0--x, 1--y, 2--z \n");            // 0
            
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
    numparsT = 2 * xyzdim[0] * xyzdim[1] + 2 * xyzdim[0] * (xyzdim[2]-2)
             + 2 * (xyzdim[1]-2) * (xyzdim[2]-2);
    
    xyzdd[0] = (xyzminmax[1] - xyzminmax[0]) / (xyzdim[0] - 1);
    xyzdd[1] = (xyzminmax[3] - xyzminmax[2]) / (xyzdim[1] - 1);
    xyzdd[2] = (xyzminmax[5] - xyzminmax[4]) / (xyzdim[2] - 1);
    

    numparsTloc = numparsT;
    numparsSloc = numparsS;
    
    memcpy(xyzdimloc, xyzdim, sizeof(xyzdim));
    memcpy(xyzminmaxloc, xyzminmax, sizeof(xyzminmax));
    memcpy(xyzddloc, xyzdd, sizeof(xyzdd));
    

    if (sflag == 2) {

        numparsSloc = numparsS / p;
        if (rank == p-1) numparsSloc += (numparsS - (numparsS / p) * p);
        
        make_vector(xS,numparsSloc);
        make_vector(yS,numparsSloc);
        make_vector(zS,numparsSloc);
        make_vector(qS,numparsSloc);
        make_vector(tenergy,numparsT);
    
        for (i = 0; i < numparsT; i++) tenergy[i] = 0.0;
        
    /* Reading in coordinates and charges for the source particles*/
        MPI_File_open(MPI_COMM_WORLD, sampin1, MPI_MODE_RDONLY, MPI_INFO_NULL, &fpmpi);
        MPI_File_seek(fpmpi, (MPI_Offset)(rank*(numparsS/p)*4*sizeof(double)), MPI_SEEK_SET);
        for (i = 0; i < numparsSloc; i++) {
            MPI_File_read(fpmpi, buf, 4, MPI_DOUBLE, &status);
            xS[i] = buf[0];
            yS[i] = buf[1];
            zS[i] = buf[2];
            qS[i] = buf[3];
        }
        MPI_File_close(&fpmpi);
        
        if (rank == 0) {
            make_vector(tenergyglob,numparsT);
            make_vector(denergyglob,numparsT);
            
            MPI_File_open(MPI_COMM_SELF, sampin3, MPI_MODE_RDONLY, MPI_INFO_NULL, &fpmpi);
            MPI_File_seek(fpmpi, (MPI_Offset)0, MPI_SEEK_SET);
            MPI_File_read(fpmpi, &time_direct, 1, MPI_DOUBLE, &status);
            MPI_File_read(fpmpi, denergyglob, numparsT, MPI_DOUBLE, &status);
            MPI_File_close(&fpmpi);
        }

        xxminmax[0][0] = xyzminmax[0];  xxminmax[0][1] = xyzminmax[1];
        xxminmax[0][2] = xyzminmax[2];  xxminmax[0][3] = xyzminmax[3];
        xxdim[0][0] = xyzdim[0];  xxdim[0][1] = xyzdim[1];
    
        xxminmax[1][0] = xyzminmax[0];  xxminmax[1][1] = xyzminmax[1];
        xxminmax[1][2] = xyzminmax[4]+xyzdd[2];  xxminmax[1][3] = xyzminmax[5]-xyzdd[2];
        xxdim[1][0] = xyzdim[0];  xxdim[1][1] = xyzdim[2]-2;
    
        xxminmax[2][0] = xyzminmax[2]+xyzdd[1];  xxminmax[2][1] = xyzminmax[3]-xyzdd[1];
        xxminmax[2][2] = xyzminmax[4]+xyzdd[2];  xxminmax[2][3] = xyzminmax[5]-xyzdd[2];
        xxdim[2][0] = xyzdim[1]-2;  xxdim[2][1] = xyzdim[2]-2;

    } else if (sflag == 0) {

        numparsTloc[0] = numparsS / p;
        if (rank == p-1) {
            numparsSloc += (numparsS - (numparsS / p) * p);
        }
        
        make_vector(xS,numparsS);
        make_vector(yS,numparsS);
        make_vector(zS,numparsS);
        make_vector(qS,numparsS);
        make_vector(tenergy,numparsT);
    
        for (i = 0; i < numparsT; i++) tenergy[i] = 0.0;
        
    /* Reading in coordinates and charges for the source particles*/
        MPI_File_open(MPI_COMM_WORLD, sampin1, MPI_MODE_RDONLY, MPI_INFO_NULL, &fpmpi);
        MPI_File_seek(fpmpi, (MPI_Offset)(rank*(numparsS/p)*4*sizeof(double)), MPI_SEEK_SET);
        for (i = 0; i < numparsSloc; i++) {
            MPI_File_read(fpmpi, buf, 4, MPI_DOUBLE, &status);
            xS[i] = buf[0];
            yS[i] = buf[1];
            zS[i] = buf[2];
            qS[i] = buf[3];
        }
        MPI_File_close(&fpmpi);
        
        if (rank == 0) {
            make_vector(tenergyglob,numparsT);
            make_vector(denergyglob,numparsT);
            
            MPI_File_open(MPI_COMM_SELF, sampin3, MPI_MODE_RDONLY, MPI_INFO_NULL, &fpmpi);
            MPI_File_seek(fpmpi, (MPI_Offset)0, MPI_SEEK_SET);
            MPI_File_read(fpmpi, &time_direct, 1, MPI_DOUBLE, &status);
            MPI_File_read(fpmpi, denergyglob, numparsT, MPI_DOUBLE, &status);
            MPI_File_close(&fpmpi);
        }

        xxminmax[0][0] = xyzminmax[0];  xxminmax[0][1] = xyzminmax[1];
        xxminmax[0][2] = xyzminmax[2];  xxminmax[0][3] = xyzminmax[3];
        xxdim[0][0] = xyzdim[0];  xxdim[0][1] = xyzdim[1];
    
        xxminmax[1][0] = xyzminmax[0];  xxminmax[1][1] = xyzminmax[1];
        xxminmax[1][2] = xyzminmax[4]+xyzdd[2];  xxminmax[1][3] = xyzminmax[5]-xyzdd[2];
        xxdim[1][0] = xyzdim[0];  xxdim[1][1] = xyzdim[2]-2;
    
        xxminmax[2][0] = xyzminmax[2]+xyzdd[1];  xxminmax[2][1] = xyzminmax[3]-xyzdd[1];
        xxminmax[2][2] = xyzminmax[4]+xyzdd[2];  xxminmax[2][3] = xyzminmax[5]-xyzdd[2];
        xxdim[2][0] = xyzdim[1]-2;  xxdim[2][1] = xyzdim[2]-2;


    }

    time2 = MPI_Wtime();
    time_preproc = time2 - time1;

    
    /* Calling main treecode subroutine to calculate approximate energy */
    treecode_grid_bdry(xS, yS, zS, qS, xxminmax[0], xxdim[0], xyzminmax[4], 0, numparsSloc,
                       tenergy, &tpengg[0], order, theta, maxparnode, time_tree,
                       pot_type, kappa);
    
    treecode_grid_bdry(xS, yS, zS, qS, xxminmax[0], xxdim[0], xyzminmax[5], 0, numparsSloc,
                       &tenergy[xxdim[0][0]*xxdim[0][1]],
                       &tpengg[1], order, theta, maxparnode, time_tree,
                       pot_type, kappa);
    
    treecode_grid_bdry(xS, yS, zS, qS, xxminmax[1], xxdim[1], xyzminmax[2], 1, numparsSloc,
                       &tenergy[2*xxdim[0][0]*xxdim[0][1]],
                       &tpengg[2], order, theta, maxparnode, time_tree,
                       pot_type, kappa);
    
    treecode_grid_bdry(xS, yS, zS, qS, xxminmax[1], xxdim[1], xyzminmax[3], 1, numparsSloc,
                       &tenergy[2*xxdim[0][0]*xxdim[0][1] + xxdim[1][0]*xxdim[1][1]],
                       &tpengg[3], order, theta, maxparnode, time_tree,
                       pot_type, kappa);
    
    treecode_grid_bdry(xS, yS, zS, qS, xxminmax[2], xxdim[2], xyzminmax[0], 2, numparsSloc,
                       &tenergy[2*xxdim[0][0]*xxdim[0][1] + 2*xxdim[1][0]*xxdim[1][1]],
                       &tpengg[4], order, theta, maxparnode, time_tree,
                       pot_type, kappa);
    
    treecode_grid_bdry(xS, yS, zS, qS, xxminmax[2], xxdim[2], xyzminmax[1], 2, numparsSloc,
                       &tenergy[2*xxdim[0][0]*xxdim[0][1] + 2*xxdim[1][0]*xxdim[1][1]
                                + xxdim[2][0]*xxdim[2][1]],
                       &tpengg[5], order, theta, maxparnode, time_tree,
                       pot_type, kappa);
    
    tpeng = sum(tpengg, 6);

    
    /* Reducing values to root process */
    MPI_Reduce(time_tree, &time_tree_glob[0], 4, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(time_tree, &time_tree_glob[1], 4, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(time_tree, &time_tree_glob[2], 4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if (rank == 0) dpengglob = sum(denergyglob, numparsT);
    MPI_Reduce(&tpeng, &tpengglob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    
    if (rank == 0) {
        /* Printing direct and treecode time calculations: */
        printf("                   Direct time (s):  %f\n\n", time_direct);
        printf("              Number of processors:  %d\n", p);
        printf("              Pre-process time (s):  %f\n", time_preproc);
        printf("      Min, Max tree setup time (s):  %f, %f\n", time_tree_glob[0][0], time_tree_glob[1][0]);
        printf("             Min, Max cp1 time (s):  %f, %f\n", time_tree_glob[0][1], time_tree_glob[1][1]);
        printf("             Min, Max cp2 time (s):  %f, %f\n", time_tree_glob[0][2], time_tree_glob[1][2]);
        
        printf("      Min, Max total tree time (s):  %f, %f\n", time_tree_glob[0][3], time_tree_glob[1][3]);
        printf(" Preproc + Max total tree time (s):  %f\n\n", time_tree_glob[1][3] + time_preproc);
    
        /* Printing error in potential energy and potential energy */
        printf("           Direct potential energy:  %f\n", dpengglob);
        printf("             Tree potential energy:  %f\n\n", tpengglob);
    
        printf("Absolute error for total potential:  %e\n",
               fabs(tpengglob-dpengglob));
        printf("Relative error for total potential:  %e\n\n",
               fabs((tpengglob-dpengglob)/dpengglob));
    }
    
    
    /* Computing pointwise potential errors */
    MPI_Reduce(tenergy, tenergyglob, numparsT, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    
    if (rank == 0) {
        
        inferr = 0.0;
        relinferr = 0.0;
        n2err = 0.0;
        reln2err = 0.0;

        for (j = 0; j < numparsT; j++) {
            temp = fabs(denergyglob[j] - tenergyglob[j]);
    
            if (temp >= inferr) inferr = temp;
            if (fabs(denergyglob[j]) >= relinferr) relinferr = fabs(denergyglob[j]);

            n2err = n2err + pow(denergyglob[j] - tenergyglob[j], 2.0);
            reln2err = reln2err + pow(denergyglob[j], 2.0);
        }

        relinferr = inferr / relinferr;
        reln2err = sqrt(n2err / reln2err);
        n2err = sqrt(n2err);
    
    
        printf("Absolute inf norm error in potential:  %e \n", inferr);
        printf("Relative inf norm error in potential:  %e \n\n", relinferr);
        printf("  Absolute 2 norm error in potential:  %e \n", n2err);
        printf("  Relative 2 norm error in potential:  %e \n\n", reln2err);
    
        fp = fopen(sampout, "a");
        fprintf(fp, "%s \t %s \t %d \t %f \t %f \t %f \t"
                    "%f \t %f \t %f \t %d \t %d \t %d \t"
                    "%f \t %d \t %d \t %f \t %d \t %d \t %d \t %d \t"
                    "%f \t %f \t %f \t %f \t %f \t %f \t %f \t"
                    "%f \t %f \t %f \t %f \t %f \t %f \t %f \t"
                    "%e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \n",
                sampin1, sampin3, numparsS, xyzminmax[0], xyzminmax[1], xyzminmax[2], //1 ends
                xyzminmax[3], xyzminmax[4], xyzminmax[5], xyzdim[0], xyzdim[1], xyzdim[2], //2 ends
                theta, order, maxparnode, kappa, pot_type, sflag, dflag, p, //3 ends
                time_preproc,
                time_tree_glob[0][0], time_tree_glob[1][0],
                time_tree_glob[2][0]/(double)p,
                time_tree_glob[0][1], time_tree_glob[1][1],
                time_tree_glob[2][1]/(double)p, //4 ends
                time_tree_glob[0][2], time_tree_glob[1][2],
                time_tree_glob[2][2]/(double)p,
                time_tree_glob[0][3], time_tree_glob[1][3],
                time_tree_glob[2][3]/(double)p,
                time_tree_glob[1][3] + time_preproc, //5 ends
                dpengglob, tpengglob, fabs(tpengglob-dpengglob),
                fabs((tpengglob-dpengglob)/dpengglob),
                inferr, relinferr, n2err, reln2err); //6 ends
        fclose(fp);
    }
    
    
    free_vector(xS);
    free_vector(yS);
    free_vector(zS);
    free_vector(qS);
    
    free_vector(denergyglob);
    free_vector(tenergyglob);
    free_vector(tenergy);
    
    MPI_Finalize();
    return 0;
    
}
