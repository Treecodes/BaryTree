#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>

#include "array.h"
#include "tools.h"


void direct_eng(double *xS, double *yS, double *zS, double *qS, 
                double *xT, double *yT, double *zT,
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
    double *xS, *yS, *zS, *qS;  /* source particles */
    double *xT, *yT, *zT;       /* target particles */
    double *denergy;  /* exact energy, treecode energy */

    /* for potential energy calculation */
    double dpeng, dpengglob;

    double kappa;

    // insert variables for date-time calculation?
    double time1, time2, time_direct, time_direct_tot;
    double time_direct_max, time_direct_min;
    
    /* input and output files */
    char *sampin1, *sampin2, *sampout;

    //local variables
    int i;
    int numparsTloc, maxparsTloc, globparsTloc;
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
            printf("       sampin2:  targets input file \n");                // "T1000000.bin"
            printf("       sampout:  direct calc potential output file \n"); // "ex_s4_t6.bin"
            printf("      numparsS:  number of sources \n");                 // 10000
            printf("      numparsT:  number of targets \n");                 // 1000000
            printf("         kappa:  screened Coulomb parameter \n");        // 0.00
            printf("      pot_type:  0--Coulomb, 1--screened Coulomb \n");   // 1
        }
        return 0;
    }
    
    sampin2 = argv[2];
    sampout = argv[3];
    numparsS = atoi(argv[4]);
    numparsT = atoi(argv[5]);
    kappa = atof(argv[6]);
    pot_type = atoi(argv[7]);
    
    
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


    /* Calling main treecode subroutine to calculate approximate energy */

    time1 = MPI_Wtime();

    direct_eng(xS, yS, zS, qS, xT, yT, zT, numparsS, numparsTloc,
                denergy, &dpeng, pot_type, kappa);

    time2 = MPI_Wtime();
    time_direct = time2-time1;

    
    /* Reducing values to root process */
    MPI_Reduce(&time_direct, &time_direct_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&time_direct, &time_direct_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&time_direct, &time_direct_tot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&dpeng, &dpengglob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    

    MPI_File_open(MPI_COMM_WORLD, sampout, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fpmpi);
    if (rank == 0)
    {
        MPI_File_seek(fpmpi, (MPI_Offset)0, MPI_SEEK_SET);
        MPI_File_write(fpmpi, &time_direct_max, 1, MPI_DOUBLE, &status);
    }
    MPI_File_seek(fpmpi, (MPI_Offset)(globparsTloc+1)*sizeof(double), MPI_SEEK_SET);
    MPI_File_write(fpmpi, denergy, numparsTloc, MPI_DOUBLE, &status);
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
    
    MPI_Finalize();
    return 0;
}


void direct_eng(double *xS, double *yS, double *zS, double *qS, 
                double *xT, double *yT, double *zT,
                int numparsS, int numparsT, double *denergy, double *dpeng,
                int pot_type, double kappa)
{
        /* local variables */
        int i, j;
        double tx, ty, tz, xi, yi, zi, teng, rad;

        *dpeng = 0.0;

        if (pot_type == 0) {
                for (i = 0; i < numparsT; i++) {
                        xi = xT[i];
                        yi = yT[i];
                        zi = zT[i];
                        teng = 0.0;

                        for (j = 0; j < numparsS; j++) {
                                tx = xi - xS[j];
                                ty = yi - yS[j];
                                tz = zi - zS[j];
                                teng = teng + qS[j] / sqrt(tx*tx + ty*ty + tz*tz);
                        }
                        denergy[i] = teng;
                }

        } else if (pot_type == 1) {
                for (i = 0; i < numparsT; i++) {
                        xi = xT[i];
                        yi = yT[i];
                        zi = zT[i];
                        teng = 0.0;

                        for (j = 0; j < numparsS; j++) {
                                tx = xi - xS[j];
                                ty = yi - yS[j];
                                tz = zi - zS[j];
                                rad = sqrt(tx*tx + ty*ty + tz*tz);
                                teng = teng + qS[j] * exp(-kappa * rad) / rad;
                        }
                        denergy[i] = teng;
                }
        }


        *dpeng = sum(denergy, numparsT);


        return;

}
