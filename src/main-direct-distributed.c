#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <omp.h>
#include <float.h>


#include "array.h"
#include "tools.h"


void direct_eng(double *xS, double *yS, double *zS, double *qS, double *wS,
                double *xT, double *yT, double *zT, double *qT,
                int numparsS, int numparsT, double *denergy, double *dpeng,
                int pot_type, double kappa, int numDevices, int numThreads);

/* The treedriver routine in Fortran */
int main(int argc, char **argv)
{

    int rank, numProcs;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    printf("rank = %i\n",rank);
    /* runtime parameters */
    int numparsS, numparsT, numDevices, numThreads;
    int pot_type;

    /* arrays for coordinates, charges, energy of target particles */
    /* source particles */
    double *xS = NULL; 
    double *yS = NULL;
    double *zS = NULL;
    double *qS = NULL;
    double *wS = NULL;

    /* target particles */
    double *xT = NULL;
    double *yT = NULL;
    double *zT = NULL;
    double *qT = NULL;

    /* exact energy */
    double *denergy = NULL;

    /* for potential energy calculation */
    double dpeng, dpengglob;

    double kappa;

    // insert variables for date-time calculation?
    double time1, time2, time_direct, time_direct_tot, timeCommunicate, timeInteract;
    double time_direct_max, time_direct_min;
    
    /* input and output files */
    char *sampin1, *sampin2, *sampout, *sampdatout;
    FILE *fp;

    //local variables
    int i;
    int numparsTloc, maxparsTloc, globparsTloc;
    int numparsSloc, maxparsSloc, globparsSloc;
    double buf[5];
    
    /* MPI Variables */
    MPI_File fpmpi;
    MPI_Status status;
    if (rank==0) printf("Setup MPI file variables.\n");
//    omp_set_num_threads(4);


    sampin1 = argv[1];
    if (strcmp(sampin1,"--help") == 0)
    {
        if (rank == 0)
        {
            printf("Input arguments: \n");
            printf("       sampin1:  sources input file \n");                // "S10000.bin"
            printf("       sampin2:  targets input file \n");                // "T1000000.bin"
            printf("       sampout:  direct calc potential output file \n"); // "ex_s4_t6.bin"
            printf("    sampdatout:  human readable data output file \n");   // "out.txt"
            printf("      numparsS:  number of sources \n");                 // 10000
            printf("      numparsT:  number of targets \n");                 // 1000000
            printf("         kappa:  screened Coulomb parameter \n");        // 0.00
            printf("      pot_type:  0--Coulomb, 1--screened Coulomb \n");   // 1
            printf("      number of devices: \n");   // 1

        }
        return 0;
    }
    
    sampin2 = argv[2];
    sampout = argv[3];
    sampdatout = argv[4];
    numparsS = atoi(argv[5]);
    numparsT = atoi(argv[6]);
    kappa = atof(argv[7]);
    pot_type = atoi(argv[8]);
    numDevices = atoi(argv[9]);
    numThreads = atoi(argv[10]);
    
    numparsTloc = (int)floor((double)numparsT/(double)numProcs);
    maxparsTloc = numparsTloc + (numparsT - (int)floor((double)numparsT/(double)numProcs) * numProcs);

	numparsSloc = (int)floor((double)numparsS/(double)numProcs);
	maxparsSloc = numparsSloc + (numparsS - (int)floor((double)numparsS/(double)numProcs) * numProcs);
    
    if (rank == 0){
    	numparsTloc = maxparsTloc;
    	numparsSloc = maxparsSloc;
    }
    globparsTloc = maxparsTloc + numparsTloc * (rank-1);
    globparsSloc = maxparsSloc + numparsSloc * (rank-1);

    printf("Proc %i has %i particles.\n", rank, numparsTloc);
    if (rank==0) printf("Max local T and local S: %i, %i.\n", maxparsTloc, maxparsSloc);

//    printf("globalparsTloc = %i\n", globparsTloc);
//    printf("globalparsSloc = %i\n", globparsSloc);


    make_vector(xS,numparsSloc);
    make_vector(yS,numparsSloc);
    make_vector(zS,numparsSloc);
    make_vector(qS,numparsSloc);
    make_vector(wS,numparsSloc);

    make_vector(xT,numparsTloc);
    make_vector(yT,numparsTloc);
    make_vector(zT,numparsTloc);
    make_vector(qT,numparsTloc);

    make_vector(denergy,numparsTloc);

    // Allocate and zero out receive buffers.
	int sendTo, recvFrom;
	double * xS_foreign;
	double * yS_foreign;
	double * zS_foreign;
	double * qS_foreign;
	double * wS_foreign;
	make_vector(xS_foreign,maxparsSloc);
	make_vector(yS_foreign,maxparsSloc);
	make_vector(zS_foreign,maxparsSloc);
	make_vector(qS_foreign,maxparsSloc);
	make_vector(wS_foreign,maxparsSloc);


    /* Reading in coordinates and charges for the source particles*/
    MPI_File_open(MPI_COMM_WORLD, sampin1, MPI_MODE_RDONLY, MPI_INFO_NULL, &fpmpi);
    MPI_File_seek(fpmpi, (MPI_Offset)globparsSloc*5*sizeof(double), MPI_SEEK_SET);
    for (i = 0; i < numparsSloc; i++) {
        MPI_File_read(fpmpi, buf, 5, MPI_DOUBLE, &status);
        xS[i] = buf[0];
        yS[i] = buf[1];
        zS[i] = buf[2];
        qS[i] = buf[3];
        wS[i] = buf[4];
    }
    MPI_File_close(&fpmpi);
    

    /* Reading in coordinates for the targets */
    MPI_File_open(MPI_COMM_WORLD, sampin2, MPI_MODE_RDONLY, MPI_INFO_NULL, &fpmpi);
    MPI_File_seek(fpmpi, (MPI_Offset)globparsTloc*4*sizeof(double), MPI_SEEK_SET);
    for (i = 0; i < numparsTloc; i++) {
        MPI_File_read(fpmpi, buf, 4, MPI_DOUBLE, &status);
        xT[i] = buf[0];
        yT[i] = buf[1];
        zT[i] = buf[2];
        qT[i] = buf[3];
    }
    MPI_File_close(&fpmpi);

    dpeng=0.0;
    for (i=0;i<maxparsTloc;i++){
		denergy[i]=0.0;
	}


    /* Interact with self */
    time1 = MPI_Wtime();
//    direct_eng(xS, yS, zS, qS, wS, xT, yT, zT, qT, numparsSloc, numparsTloc,
//                denergy, &dpeng, pot_type, kappa, numDevices, numThreads);

//    if (rank==0) printf("Time to compute self interaction: %f\n", MPI_Wtime()-time1);
//    printf("Proc %i: Time to compute self interaction: %f\n", rank, MPI_Wtime()-time1);

    /* Get data from other processors, interact with them. */



	MPI_Request request1s, request2s, request3s, request4s, request5s;
	MPI_Request request1r, request2r, request3r, request4r, request5r;

	for (i=0;i<maxparsSloc;i++){
		xS_foreign[i]=0.0;
		yS_foreign[i]=0.0;
		zS_foreign[i]=0.0;
		qS_foreign[i]=0.0;
		wS_foreign[i]=0.0;
	}
    for (int procID=1;procID<numProcs;procID++){
    	timeCommunicate = MPI_Wtime();
//    	MPI_Barrier(MPI_COMM_WORLD);

    	// Send and receive source particles
    	sendTo = (rank+procID)%numProcs;
    	recvFrom = (numProcs+rank-procID)%numProcs;

    	MPI_Isend(xS, numparsSloc, MPI_DOUBLE, sendTo, 1, MPI_COMM_WORLD, &request1s);
    	MPI_Isend(yS, numparsSloc, MPI_DOUBLE, sendTo, 2, MPI_COMM_WORLD, &request2s);
    	MPI_Isend(zS, numparsSloc, MPI_DOUBLE, sendTo, 3, MPI_COMM_WORLD, &request3s);
    	MPI_Isend(qS, numparsSloc, MPI_DOUBLE, sendTo, 4, MPI_COMM_WORLD, &request4s);
    	MPI_Isend(wS, numparsSloc, MPI_DOUBLE, sendTo, 5, MPI_COMM_WORLD, &request5s);

    	MPI_Irecv(xS_foreign, maxparsSloc, MPI_DOUBLE, recvFrom, 1, MPI_COMM_WORLD, &request1r);
    	MPI_Irecv(yS_foreign, maxparsSloc, MPI_DOUBLE, recvFrom, 2, MPI_COMM_WORLD, &request2r);
    	MPI_Irecv(zS_foreign, maxparsSloc, MPI_DOUBLE, recvFrom, 3, MPI_COMM_WORLD, &request3r);
    	MPI_Irecv(qS_foreign, maxparsSloc, MPI_DOUBLE, recvFrom, 4, MPI_COMM_WORLD, &request4r);
    	MPI_Irecv(wS_foreign, maxparsSloc, MPI_DOUBLE, recvFrom, 5, MPI_COMM_WORLD, &request5r);

//    	if (procID==1){
//    		// Overlap communication and computation.
//    		direct_eng(xS, yS, zS, qS, wS, xT, yT, zT, qT, numparsSloc, numparsTloc,
//    		                denergy, &dpeng, pot_type, kappa, numDevices, numThreads);
//    	}


    	direct_eng(xS, yS, zS, qS, wS, xT, yT, zT, qT, maxparsSloc, numparsTloc,
    					denergy, &dpeng, pot_type, kappa, numDevices, numThreads);

    	MPI_Wait(&request1r, &status);
    	MPI_Wait(&request2r, &status);
    	MPI_Wait(&request3r, &status);
    	MPI_Wait(&request4r, &status);
    	MPI_Wait(&request5r, &status);

//    	direct_eng(xS, yS, zS, qS, wS, xT, yT, zT, qT, maxparsSloc, numparsTloc,
//				denergy, &dpeng, pot_type, kappa, numDevices, numThreads);

    	MPI_Wait(&request1s, &status);
		MPI_Wait(&request2s, &status);
		MPI_Wait(&request3s, &status);
		MPI_Wait(&request4s, &status);
		MPI_Wait(&request5s, &status);


//		direct_eng(xS, yS, zS, qS, wS, xT, yT, zT, qT, maxparsSloc, numparsTloc,
//				denergy, &dpeng, pot_type, kappa, numDevices, numThreads);

//    	if (rank==0) printf("Time to communicate: %f\n", MPI_Wtime()-timeCommunicate);
//    	printf("Proc %i: Time to communicate: %f\n", rank, MPI_Wtime()-timeCommunicate);


//    	MPI_Barrier(MPI_COMM_WORLD);

//    	printf("Proc %i: Send to %i, Recv from %i\n", rank, sendTo, recvFrom);

    	// Compute interaction
//    	timeInteract = MPI_Wtime();
//    	direct_eng(xS_foreign, yS_foreign, zS_foreign, qS_foreign, wS_foreign, xT, yT, zT, qT, maxparsSloc, numparsTloc,
//				denergy, &dpeng, pot_type, kappa, numDevices, numThreads);
//    	printf("Proc %i: Time to compute interaction: %f\n", rank, MPI_Wtime()-timeInteract);


//    	MPI_Barrier(MPI_COMM_WORLD);

    	for (i=0;i<maxparsSloc;i++){ // re-zero out the buffers.  Maybe necessary since one proc has more data than rest?
			xS[i] = xS_foreign[i];
			yS[i] = yS_foreign[i];
			zS[i] = zS_foreign[i];
			qS[i] = qS_foreign[i];
			wS[i] = wS_foreign[i];
			if (procID<numProcs-1){
				xS_foreign[i]=0.0;
				yS_foreign[i]=0.0;
				zS_foreign[i]=0.0;
				qS_foreign[i]=0.0;
				wS_foreign[i]=0.0;
			}
		}

//    	printf("Proc %i, dpeng = %f\n",rank,dpeng);
    }

    direct_eng(xS, yS, zS, qS, wS, xT, yT, zT, qT, maxparsSloc, numparsTloc,
    				denergy, &dpeng, pot_type, kappa, numDevices, numThreads);

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
        printf("                 Avg direct time (s):  %f\n", time_direct_tot/(double)numProcs);
        printf("   Total direct time (s) on %d procs:  %f\n\n", numProcs, time_direct_tot);
    
        /* Calculating value dpeng by summing all values in denergy */
        printf("             Direct potential energy:  %.15f\n\n\n", dpengglob);
    }

    
    if (rank == 0) {
        fp = fopen(sampdatout, "a");
        fprintf(fp, "%s \t %s \t %s \t %d \t %d \t %f \t %d \t %d \t"
                "%f \t %f \t %f \t %e \n",
                sampin1, sampin2, sampout, numparsS, numparsT,
                kappa, pot_type, numProcs, time_direct_max, time_direct_min,
                time_direct_tot/(double)numProcs, dpengglob);
        fclose(fp);
    }

    free_vector(xS);
    free_vector(yS);
    free_vector(zS);
    free_vector(qS);
    free_vector(wS);
    free_vector(xT);
    free_vector(yT);
    free_vector(zT);
    free_vector(denergy);
    free_vector(xS_foreign);
	free_vector(yS_foreign);
	free_vector(zS_foreign);
	free_vector(qS_foreign);
	free_vector(wS_foreign);

    MPI_Finalize();
    return 0;
}


void direct_eng( double *xS, double *yS, double *zS, double *qS, double *wS,
				double *xT, double *yT, double *zT, double *qT,
                int numparsS, int numparsT, double *denergy, double *dpeng,
                int pot_type, double kappa, int numDevices, int numThreads)
{
        /* local variables */
        int i, j;
        double tx, ty, tz, xi, yi, zi, qi, teng, rad;
//        printf("numparsT,  numparsS: %i, %i\n", numparsT, numparsS);

		if (pot_type == 0 || pot_type == 4) { // Coulomb with singularity skipping.  Lagrange or Hermite.


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
						if (rad>1e-14){
							teng = teng + qS[j]*wS[j] / rad;
						}
				}
				denergy[i] +=  teng;

			}


//        	queue = (queue%2)+1;


		} else if (pot_type == 1 || pot_type == 5) {
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
								if (rad>1e-14){
									teng = teng + qS[j]*wS[j] * exp(-kappa * rad) / rad;
								}
						}
						denergy[i] += teng;
				}

		} else if (pot_type == 2 || pot_type == 6) {
			double kappaSq = kappa*kappa;
				for (i = 0; i < numparsT; i++) {
						xi = xT[i];
						yi = yT[i];
						zi = zT[i];
						qi = qT[i];
						teng = 2*M_PI*kappaSq*qi;  // 2pi alpha^2*f_t for SS scheme exp(-r^2/alpha^2)

						for (j = 0; j < numparsS; j++) {
								tx = xi - xS[j];
								ty = yi - yS[j];
								tz = zi - zS[j];
								rad = sqrt(tx*tx + ty*ty + tz*tz);
								if (rad>1e-14){
									teng = teng + ( qS[j] - qi* exp(-rad*rad/kappaSq)) *wS[j]/ rad;
								}
						}
						denergy[i] += teng;
				}

		} else if (pot_type == 3 || pot_type == 7) {

			for (i = 0; i < numparsT; i++) {
					xi = xT[i];
					yi = yT[i];
					zi = zT[i];
					qi = qT[i];
					teng = 4*M_PI*qi/kappa/kappa;  // 4pi*f_t/k^2

					for (j = 0; j < numparsS; j++) {
							tx = xi - xS[j];
							ty = yi - yS[j];
							tz = zi - zS[j];
							rad = sqrt(tx*tx + ty*ty + tz*tz);
							if (rad>1e-14){
								teng += ( qS[j] - qi) *wS[j] * exp(-kappa * rad) / rad;
							}
					}
					denergy[i] += teng;
			}

		} // end pot=3 or 7




        *dpeng = sum(denergy, numparsT);


//} // end pragma omp parallel


        return;

}
