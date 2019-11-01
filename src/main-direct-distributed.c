#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <float.h>


#include "array.h"
#include "tools.h"
#include "direct.h"
#include "kernels/kernels.h"

int main(int argc, char **argv)
{
    int rank, numProcs, provided;
    
    MPI_Init(&argc, &argv);
    int ierr;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);



    /* runtime parameters */
    int numparsS, numparsT;
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
    double time1, time2, time3, time4, time_direct, time_direct_tot, timeCommunicate, timeInteract;
    double time_direct_max, time_direct_min;
    double total_time_start, total_time_stop;
    
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

    sampin1 = argv[1];
    if (strcmp(sampin1,"--help") == 0)
    {
        if (rank == 0)
        {
            printf("Input arguments: \n");
            printf("       infile 1:  sources input file \n");                // "S10000.bin"
            printf("       infile 2:  targets input file \n");                // "T1000000.bin"
            printf("    pot outfile:  direct calc potential output file \n"); // "ex_s4_t6.bin"
            printf("    csv outfile:  human readable data output file \n");   // "out.txt"
            printf("       numparsS:  number of sources \n");                 // 10000
            printf("       numparsT:  number of targets \n");                 // 1000000
            printf("       pot type:  0--Coulomb\n");
            printf("                  1--screened Coulomb/Yukawa\n");
            printf("          kappa:  screened Coulomb parameter \n");        // 0.00


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
    

// Set up kernel
    double (*kernel)( double targetX, double targetY, double targetZ, double targetQ,
    			double sourceX, double sourceY, double sourceZ, double sourceQ, double sourceW,
				double kappa);

    if (pot_type==0){
    	kernel = &coulombKernel;
    }else{
    	if (pot_type==1){
    		kernel = &yukawaKernel;
    	}

    }

    numparsTloc = (int)floor((double)numparsT/(double)numProcs);
    maxparsTloc = numparsTloc + (numparsT - (int)floor((double)numparsT/(double)numProcs) * numProcs);

	numparsSloc = (int)floor((double)numparsS/(double)numProcs);
	maxparsSloc = numparsSloc + (numparsS - (int)floor((double)numparsS/(double)numProcs) * numProcs);
    if (rank==0) printf("Num local T and local S: %i, %i.\n", numparsSloc, numparsTloc);
    
    if (rank == 0){
    	numparsTloc = maxparsTloc;
    	numparsSloc = maxparsSloc;
    }

    globparsTloc = maxparsTloc + numparsTloc * (rank-1);
    globparsSloc = maxparsSloc + numparsSloc * (rank-1);

    if (rank==0) printf("Max local T and local S: %i, %i.\n", maxparsTloc, maxparsSloc);

    double *S_local;
    make_vector(S_local, 5*maxparsSloc);
    xS = &S_local[0*maxparsSloc];
    yS = &S_local[1*maxparsSloc];
    zS = &S_local[2*maxparsSloc];
    qS = &S_local[3*maxparsSloc];
    wS = &S_local[4*maxparsSloc];


    for (int i = 0; i < 5*maxparsSloc; i++) {
    	S_local[i] = 0.0;
    }


    double *T_local;
    make_vector(T_local, 4*numparsTloc);
    xT = &T_local[0*numparsTloc];
    yT = &T_local[1*numparsTloc];
    zT = &T_local[2*numparsTloc];
    qT = &T_local[3*numparsTloc];


    make_vector(denergy,numparsTloc);

    // Allocate and zero out receive buffers.
	int sendTo, recvFrom;
	double *xS_foreign1;
	double *yS_foreign1;
	double *zS_foreign1;
	double *qS_foreign1;
	double *wS_foreign1;
	double *S_foreign1;
	double *xS_foreign2;
	double *yS_foreign2;
	double *zS_foreign2;
	double *qS_foreign2;
	double *wS_foreign2;
	double *S_foreign2;

	make_vector(S_foreign1, 5*maxparsSloc);
//	memset(S_foreign1,0,5*maxparsSloc);
	xS_foreign1 = &S_foreign1[0*maxparsSloc];
	yS_foreign1 = &S_foreign1[1*maxparsSloc];
	zS_foreign1 = &S_foreign1[2*maxparsSloc];
	qS_foreign1 = &S_foreign1[3*maxparsSloc];
	wS_foreign1 = &S_foreign1[4*maxparsSloc];

	make_vector(S_foreign2, 5*maxparsSloc);
	xS_foreign2 = &S_foreign2[0*maxparsSloc];
	yS_foreign2 = &S_foreign2[1*maxparsSloc];
	zS_foreign2 = &S_foreign2[2*maxparsSloc];
	qS_foreign2 = &S_foreign2[3*maxparsSloc];
	wS_foreign2 = &S_foreign2[4*maxparsSloc];


    /* Reading in coordinates and charges for the source particles*/
    ierr = MPI_File_open(MPI_COMM_WORLD, sampin1, MPI_MODE_RDONLY, MPI_INFO_NULL, &fpmpi);
    if (ierr != MPI_SUCCESS) {
        fprintf(stderr,"FILE COULD NOT OPEN\n");
        exit(1);
    }

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
    
//    MPI_Barrier(MPI_COMM_WORLD);

    /* Reading in coordinates for the targets */
    ierr = MPI_File_open(MPI_COMM_WORLD, sampin2, MPI_MODE_RDONLY, MPI_INFO_NULL, &fpmpi);
    if (ierr != MPI_SUCCESS) {
        fprintf(stderr,"FILE COULD NOT OPEN\n");
        exit(1);
    }

    MPI_File_seek(fpmpi, (MPI_Offset)globparsTloc*4*sizeof(double), MPI_SEEK_SET);
    for (i = 0; i < numparsTloc; i++) {
        MPI_File_read(fpmpi, buf, 4, MPI_DOUBLE, &status);
        xT[i] = buf[0];
        yT[i] = buf[1];
        zT[i] = buf[2];
        qT[i] = buf[3];
    }
    MPI_File_close(&fpmpi);

    dpeng = 0.0;
    for (i = 0; i < maxparsTloc; i++) {
		denergy[i] = 0.0;
	}
 
#ifdef OPENACC_ENABLED
    #pragma acc set device_num(rank) device_type(acc_device_nvidia)
    #pragma acc init device_type(acc_device_nvidia)
#endif

    /* Interact with self */
    time1 = MPI_Wtime();
    total_time_start = MPI_Wtime();


	MPI_Request request1s;
	MPI_Request request1r;

	for (i=0; i < 5*maxparsSloc; i++) {
		S_foreign1[i]=0.0;
		S_foreign2[i]=0.0;
	}


	if (numProcs == 1) { // one 1 proc, won't enter into the round robin below.  So just interact with self here.
		direct_eng(xS, yS, zS, qS, wS, xT, yT, zT, qT, maxparsSloc, numparsTloc,
						denergy, &dpeng, pot_type, kappa, (*kernel));
	} else {

		for (int procID = 1; procID < numProcs; procID++) {
			timeCommunicate = MPI_Wtime();

			// Send and receive source particles
			sendTo = (rank+procID)%numProcs;
			recvFrom = (numProcs+rank-procID)%numProcs;

			if (procID%2==1) { // communicate w/ foreign1, compute on foreign2
				MPI_Isend(S_local, 5*maxparsSloc, MPI_DOUBLE, sendTo, 1, MPI_COMM_WORLD, &request1s);
				MPI_Irecv(S_foreign1, 5*maxparsSloc, MPI_DOUBLE, recvFrom, 1, MPI_COMM_WORLD, &request1r);


				if (procID==1) { // in first iteration of loop, interact with self.
					direct_eng(xS, yS, zS, qS, wS, xT, yT, zT, qT, maxparsSloc, numparsTloc,
											denergy, &dpeng, pot_type, kappa, (*kernel));
				} else { // in subsequent iterations, interact with foreign
					direct_eng(xS_foreign2, yS_foreign2, zS_foreign2, qS_foreign2, wS_foreign2, xT, yT, zT, qT, maxparsSloc, numparsTloc, denergy, &dpeng, pot_type, kappa, (*kernel));
					for (i=0;i<5*maxparsSloc;i++){
							S_foreign2[i]=0.0;
						}
				}
			} else { // procId%2=0, communicate foreign2 and compute with foreign1
				MPI_Isend(S_local, 5*maxparsSloc, MPI_DOUBLE, sendTo, 1, MPI_COMM_WORLD, &request1s);
				MPI_Irecv(S_foreign2, 5*maxparsSloc, MPI_DOUBLE, recvFrom, 1, MPI_COMM_WORLD, &request1r);

				direct_eng(xS_foreign1, yS_foreign1, zS_foreign1, qS_foreign1, wS_foreign1, xT, yT, zT, qT, maxparsSloc, numparsTloc, denergy, &dpeng, pot_type, kappa, (*kernel));
				for (i=0;i<5*maxparsSloc;i++){
						S_foreign1[i]=0.0;
					}
			}

			time4 = MPI_Wtime();
			MPI_Wait(&request1r, &status);
			MPI_Wait(&request1s, &status);



		} // end round robin

		if ((numProcs-1)%2==1) { // in final loop, S_foreign1 was received but not yet computed with
			direct_eng(xS_foreign1, yS_foreign1, zS_foreign1, qS_foreign1, wS_foreign1, xT, yT, zT, qT, maxparsSloc, numparsTloc,
											denergy, &dpeng, pot_type, kappa, (*kernel));
		} else { // S_foreign2 is the one that needs to be computed with
			direct_eng(xS_foreign2, yS_foreign2, zS_foreign2, qS_foreign2, wS_foreign2, xT, yT, zT, qT, maxparsSloc, numparsTloc,
											denergy, &dpeng, pot_type, kappa, (*kernel));
		}
	}

    time2 = MPI_Wtime();
    time_direct = time2-time1;

    
    /* Reducing values to root process */
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&time_direct, &time_direct_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&time_direct, &time_direct_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&time_direct, &time_direct_tot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&dpeng, &dpengglob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
//        } // end OMP PARALLEL REGION
    total_time_stop = MPI_Wtime();


    MPI_File_open(MPI_COMM_WORLD, sampout, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fpmpi);
    if (rank == 0) {
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
        printf("		          Wallclock time (s): %f\n", total_time_stop-total_time_start);
        printf("		         Cumulative time (s): %f\n", (total_time_stop-total_time_start)*numProcs);
        printf("             Direct potential energy:  %.15f\n\n\n", dpengglob);
    }

    
    if (rank == 0) {
        fp = fopen(sampdatout, "a");
        fprintf(fp, "%s, %s, %s, %d, %d, %f, %d, %d,"
                "%f, %f, %f, %e \n",
                sampin1, sampin2, sampout, numparsS, numparsT,
                kappa, pot_type, numProcs, time_direct_max, time_direct_min,
                time_direct_tot/(double)numProcs, dpengglob);
        fclose(fp);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank==0) printf("Finished writing files.\n");
    free_vector(T_local);
    free_vector(S_local);


    free_vector(denergy);

	free_vector(S_foreign1);
	free_vector(S_foreign2);
    MPI_Barrier(MPI_COMM_WORLD);
	if (rank==0) printf("Finished freeing memory.\n");
    MPI_Finalize();
    return 0;
}
