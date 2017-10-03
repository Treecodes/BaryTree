#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>

#include "array.h"
#include "tools.h"


void direct_eng(double *xS, double *yS, double *zS, double *qS, 
                double *xyzminmax, int *xyzdim, double *xyzdd,
                int numparsS, int numparsT, double *denergy, double *dpeng,
                int pot_type, double kappa);

/* The treedriver routine in Fortran */
int main(int argc, char **argv)
{

    int rank, p;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    
    /* runtime parameters */
    int numparsS, numparsT;
    int pot_type;

    /* arrays for coordinates, charges, energy of target particles */
    /* source particles */
    double *xS = NULL; 
    double *yS = NULL;
    double *zS = NULL;
    double *qS = NULL;

    /* target particles */
    double xyzminmax[6];
    double xyzdd[3];
    int xyzdim[3];

    /* exact energy */
    double *denergy = NULL;

    /* for potential energy calculation */
    double dpeng, dpengglob;

    double kappa;

    // insert variables for date-time calculation?
    double time1, time2, time_direct, time_direct_tot;
    double time_direct_max, time_direct_min;
    
    /* input and output files */
    char *sampin1, *sampin2=NULL, *sampout, *sampdatout;
    FILE *fp;

    //local variables
    int i;
    int numparsSloc, maxparsSloc, globparsSloc;
    double buf[5];
    
    /* MPI Variables */
    MPI_File fpmpi;
    MPI_Status status;


    sampin1 = argv[1];
    if (strcmp(sampin1,"--help") == 0)
    {
        if (rank == 0)
        {
            printf("Input arguments: \n");
            printf("       sampin1:  sources input file \n");                // "S10000.bin"
            printf("       sampout:  direct calc potential output file \n"); // "ex_s4_t6.bin"
            printf("    sampdatout:  human readable data output file \n");   // "out.txt"
            printf("      numparsS:  number of sources \n");                 // 10000
            printf("         kappa:  screened Coulomb parameter \n");        // 0.00
            printf("      pot_type:  0--Coulomb, 1--screened Coulomb \n");   // 1

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
    
    sampin1 = argv[1];
    sampout = argv[2];
    sampdatout = argv[3];
    numparsS = atoi(argv[4]);
    kappa = atof(argv[5]);
    pot_type = atoi(argv[6]);

    xyzminmax[0] = atof(argv[7]);
    xyzminmax[1] = atof(argv[8]);
    xyzminmax[2] = atof(argv[9]);
    xyzminmax[3] = atof(argv[10]);
    xyzminmax[4] = atof(argv[11]);
    xyzminmax[5] = atof(argv[12]);

    xyzdim[0] = atoi(argv[13]);
    xyzdim[1] = atoi(argv[14]);
    xyzdim[2] = atoi(argv[25]);


    numparsT = xyzdim[0] * xyzdim[1] * xyzdim[2];
    
    xyzdd[0] = (xyzminmax[1] - xyzminmax[0]) / (xyzdim[0] - 1);
    xyzdd[1] = (xyzminmax[3] - xyzminmax[2]) / (xyzdim[1] - 1);
    xyzdd[2] = (xyzminmax[5] - xyzminmax[4]) / (xyzdim[2] - 1);
    
    
    
    numparsSloc = (int)floor((double)numparsS/(double)p);
    maxparsSloc = numparsSloc + (numparsS - (int)floor((double)numparsS/(double)p) * p);
    
    if (rank == 0) numparsSloc = maxparsSloc;
    globparsSloc = maxparsSloc + numparsSloc * (rank-1);


    make_vector(xS,numparsSloc);
    make_vector(yS,numparsSloc);
    make_vector(zS,numparsSloc);
    make_vector(qS,numparsSloc);

    make_vector(denergy,numparsT);


    /* Reading in coordinates and charges for the source particles*/
    MPI_File_open(MPI_COMM_WORLD, sampin1, MPI_MODE_RDONLY, MPI_INFO_NULL, &fpmpi);
    MPI_File_seek(fpmpi, (MPI_Offset)globparsSloc*4*sizeof(double), MPI_SEEK_SET);
    for (i = 0; i < numparsSloc; i++) {
        MPI_File_read(fpmpi, buf, 4, MPI_DOUBLE, &status);
        xS[i] = buf[0];
        yS[i] = buf[1];
        zS[i] = buf[2];
        qS[i] = buf[3];
    }
    MPI_File_close(&fpmpi);
    


    /* Calling main treecode subroutine to calculate approximate energy */

    time1 = MPI_Wtime();

    direct_eng(xS, yS, zS, qS, xyzminmax, xyzdim, xyzdd, numparsSloc, numparsT,
                denergy, &dpeng, pot_type, kappa);

    time2 = MPI_Wtime();
    time_direct = time2-time1;

    
    /* Reducing values to root process */
    MPI_Reduce(&time_direct, &time_direct_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&time_direct, &time_direct_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&time_direct, &time_direct_tot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&dpeng, &dpengglob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(denergy, denergy, numparsT, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    

    MPI_File_open(MPI_COMM_WORLD, sampout, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fpmpi);
    if (rank == 0)
    {
        MPI_File_seek(fpmpi, (MPI_Offset)0, MPI_SEEK_SET);
        MPI_File_write(fpmpi, &time_direct_max, 1, MPI_DOUBLE, &status);
        MPI_File_seek(fpmpi, (MPI_Offset)(1)*sizeof(double), MPI_SEEK_SET);
        MPI_File_write(fpmpi, denergy, numparsT, MPI_DOUBLE, &status);
    }
    MPI_File_close(&fpmpi);

    if (rank == 0)
    {
        /* Printing direct and treecode time calculations: */
        printf("                 Max direct time (s):  %f\n", time_direct_max);
        printf("                 Min direct time (s):  %f\n", time_direct_min);
        printf("                 Avg direct time (s):  %f\n", time_direct_tot/(double)p);
        printf("   Total direct time (s) on %d procs:  %f\n\n", p, time_direct_tot);
    
        /* Calculating value dpeng by summing all values in denergy */
        printf("             Direct potential energy:  %f\n", dpengglob);
    }

    
    if (rank == 0) {
        fp = fopen(sampdatout, "a");
        fprintf(fp, "%s \t %s \t %s \t %d \t %d \t %f \t %d \t %d \t"
                "%f \t %f \t %f \t %e \n",
                sampin1, sampin2, sampout, numparsS, numparsT,
                kappa, pot_type, p, time_direct_max, time_direct_min,  
                time_direct_tot/(double)p, dpengglob);
        fclose(fp);
    }

    free_vector(xS);
    free_vector(yS);
    free_vector(zS);
    free_vector(qS);
    free_vector(denergy);

    MPI_Finalize();
    return 0;
}




void direct_eng(double *xS, double *yS, double *zS, double *qS, 
                double *xyzminmax, int *xyzdim, double *xyzdd,
                int numparsS, int numparsT, double *denergy, double *dpeng,
                int pot_type, double kappa)
{
        /* local variables */
        int i, j, k, kk, nn;
        double tx, ty, tz, xi, yi, zi, teng, rad;

        nn = -1;
        *dpeng = 0.0;

        if (pot_type == 0) {
          for (i = 0; i < xyzdim[0]; i++) {
            xi = xyzminmax[0] + i * xyzdd[0];
            for (j = 0; j < xyzdim[1]; j++) {
              yi = xyzminmax[2] + j * xyzdd[1];
              for (k = 0; k < xyzdim[2]; k++) {
                zi = xyzminmax[4] + k * xyzdd[2];
                teng = 0.0;
                nn++;

                for (kk = 0; kk < numparsS; kk++) {
                    tx = xi - xS[kk];
                    ty = yi - yS[kk];
                    tz = zi - zS[kk];
                    teng = teng + qS[kk] / sqrt(tx*tx + ty*ty + tz*tz);
                }
                denergy[nn] = teng;
              }
            }
          }

        } else if (pot_type == 1) {
          for (i = 0; i < xyzdim[0]; i++) {
            xi = xyzminmax[0] + i * xyzdd[0];
            for (j = 0; j < xyzdim[1]; j++) {
              yi = xyzminmax[2] + j * xyzdd[1];
              for (k = 0; k < xyzdim[2]; k++) {
                zi = xyzminmax[4] + k * xyzdd[2];
                teng = 0.0;
                nn++;

                for (kk = 0; kk < numparsS; kk++) {
                    tx = xi - xS[kk];
                    ty = yi - yS[kk];
                    tz = zi - zS[kk];
                    rad = sqrt(tx*tx + ty*ty + tz*tz);
                    teng = teng + qS[kk] * exp(-kappa * rad) / rad;
                }
                denergy[nn] = teng;
              }
            }
          }
        }

        *dpeng = sum(denergy, numparsT);

        return;

}
